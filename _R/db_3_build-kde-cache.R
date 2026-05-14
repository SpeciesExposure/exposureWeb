# =============================================================================
# 03_build-kde-cache.R
#
# Build script – run once (locally or in CI) to pre-render per-species KDE
# temperature-distribution figures and save them as base64-encoded PNG data
# URIs in data/species-kde-cache.rds.
#
# The output is consumed by index.qmd at render time; pre-rendering here
# keeps the Quarto render fast.
#
# Usage:
#   Rscript _R/03_build-kde-cache.R
#
# Environment variables (optional):
#   DATA_DIR        – path to data directory (default: data/)
#   CACHE_DIR       – path to write output   (default: data/)
#   DEV_N_SPECIES   – integer; limit to first N species (for testing)
# =============================================================================


# =============================================================================
# 0.  Packages and helpers
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(base64enc)
  library(purrr)
})

source("_R/timeseries-utils.R")
source("_R/timeseries-plotting.R")


# =============================================================================
# 1.  Paths and config
# =============================================================================

data_dir  <- Sys.getenv("DATA_DIR",  "data")
cache_dir <- Sys.getenv("CACHE_DIR", "data")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

raster_path  <- file.path(data_dir, "temp__12.tif")
allcell_path <- file.path(data_dir, "AllCellExposureSpXVar.rds")
ts_cache_path <- file.path(data_dir, "species-timeseries-cache.rds")

for (p in c(raster_path, allcell_path, ts_cache_path)) {
  if (!file.exists(p)) stop(sprintf("Required file not found: %s", p))
}

N_GRID          <- 300L
PERIOD_BOUNDARY <- 2022L


# =============================================================================
# 2.  Load data (once, shared across all species)
# =============================================================================

cat("[kde-cache] Loading temperature raster...\n")
r <- terra::rast(raster_path)
layer_years <- as.integer(sub("WY", "", names(r)))
cat(sprintf("[kde-cache] Raster: %d layers, years %d–%d\n",
            terra::nlyr(r), min(layer_years), max(layer_years)))

cat("[kde-cache] Loading allcell data...\n")
allcell_df <- readRDS(allcell_path)
allcell_list <- split(allcell_df, allcell_df$spName)
cat(sprintf("[kde-cache] allcell: %d species\n", length(allcell_list)))
rm(allcell_df); gc()

cat("[kde-cache] Loading timeseries cache...\n")
ts_cache <- readRDS(ts_cache_path)
all_species <- names(ts_cache)

# DEV MODE
dev_n <- suppressWarnings(as.integer(Sys.getenv("DEV_N_SPECIES", "")))
if (!is.na(dev_n) && dev_n > 0) {
  all_species <- head(all_species, dev_n)
  cat(sprintf("[kde-cache] DEV MODE: limiting to first %d species\n", dev_n))
}

cat(sprintf("[kde-cache] Building KDE figures for %d species...\n", length(all_species)))


# =============================================================================
# 3.  Helper: build combined KDE + cumulative figure for one species
# =============================================================================

build_kde_figure <- function(sp_name, allcell_sp, ts_data) {

  if (is.null(allcell_sp) || nrow(allcell_sp) == 0) return(NULL)

  range_cells <- unique(allcell_sp$cell)
  range_size  <- allcell_sp$rangeSize[1]

  # --- Extract and reshape temperature values --------------------------------
  temp_mat  <- terra::extract(r, range_cells)
  temp_long <- as.data.frame(temp_mat) |>
    dplyr::mutate(cell = range_cells) |>
    tidyr::pivot_longer(cols = -cell, names_to = "layer", values_to = "temp_k") |>
    dplyr::mutate(
      year   = layer_years[match(layer, names(r))],
      temp_c = temp_k - 273.15
    ) |>
    dplyr::filter(!is.na(temp_c), year >= 1941, year <= 2025)

  if (nrow(temp_long) == 0) return(NULL)

  # --- Shared temperature grid -----------------------------------------------
  temp_range_lines <- temp_long |>
    dplyr::group_by(year) |>
    dplyr::summarise(temp_min = min(temp_c), temp_max = max(temp_c), .groups = "drop")

  temp_range <- range(temp_long$temp_c)
  grid_temps <- seq(temp_range[1], temp_range[2], length.out = N_GRID)
  grid_step  <- diff(grid_temps)[1]
  all_years  <- sort(unique(temp_long$year))

  # --- Historical ECDF and thresholds ----------------------------------------
  hist_ecdf_fn <- stats::ecdf(temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY])

  temp_binned <- tidyr::expand_grid(year = all_years, temp_mid = grid_temps) |>
    dplyr::mutate(pct_range = hist_ecdf_fn(temp_mid)) |>
    dplyr::left_join(temp_range_lines, by = "year") |>
    dplyr::filter(temp_mid >= temp_min - grid_step / 2,
                  temp_mid <= temp_max + grid_step / 2) |>
    dplyr::mutate(
      ymin = pmax(temp_mid - grid_step / 2, temp_min),
      ymax = pmin(temp_mid + grid_step / 2, temp_max)
    )

  hist_tmax <- temp_long |>
    dplyr::filter(year <= PERIOD_BOUNDARY) |>
    dplyr::group_by(cell) |>
    dplyr::summarise(cell_p99 = quantile(temp_c, 0.99), .groups = "drop") |>
    dplyr::summarise(val = quantile(cell_p99, 0.99)) |>
    dplyr::pull(val)

  hist_tmin <- temp_long |>
    dplyr::filter(year <= PERIOD_BOUNDARY) |>
    dplyr::group_by(cell) |>
    dplyr::summarise(cell_p01 = quantile(temp_c, 0.01), .groups = "drop") |>
    dplyr::summarise(val = quantile(cell_p01, 0.01)) |>
    dplyr::pull(val)

  cdf_lo    <- hist_ecdf_fn(hist_tmin)
  cdf_hi    <- hist_ecdf_fn(hist_tmax)
  hist_tmed <- median(temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY])

  # --- Individual cell lines (subsampled) ------------------------------------
  N_CELL_LINES <- 300L
  set.seed(42)
  cell_sample <- sample(range_cells, min(N_CELL_LINES, length(range_cells)))

  cell_segments <- temp_long |>
    dplyr::filter(cell %in% cell_sample) |>
    dplyr::arrange(cell, year) |>
    dplyr::group_by(cell) |>
    dplyr::mutate(
      year_end   = dplyr::lead(year),
      temp_end   = dplyr::lead(temp_c),
      cdf_val    = hist_ecdf_fn((temp_c + dplyr::coalesce(temp_end, temp_c)) / 2),
      exposed_hi = temp_c > hist_tmax,
      exposed_lo = temp_c < hist_tmin
    ) |>
    dplyr::filter(!is.na(year_end)) |>
    dplyr::ungroup()

  # --- Heatmap panel ---------------------------------------------------------
  p_kde <- ggplot(temp_binned, aes(fill = pct_range)) +
    geom_rect(aes(xmin = year - 0.5, xmax = year + 0.5, ymin = ymin, ymax = ymax)) +
    geom_line(data = temp_range_lines, aes(x = year, y = temp_min),
              inherit.aes = FALSE, colour = "black", linewidth = 1) +
    geom_line(data = temp_range_lines, aes(x = year, y = temp_max),
              inherit.aes = FALSE, colour = "black", linewidth = 1) +
    geom_segment(data = dplyr::filter(cell_segments, !exposed_hi, !exposed_lo),
                 aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
                 inherit.aes = FALSE, colour = "black", alpha = 0.12, linewidth = 0.25) +
    geom_segment(data = dplyr::filter(cell_segments, exposed_hi),
                 aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
                 inherit.aes = FALSE, colour = "#d73027", alpha = 0.6, linewidth = 0.4) +
    geom_segment(data = dplyr::filter(cell_segments, exposed_lo),
                 aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
                 inherit.aes = FALSE, colour = "#313695", alpha = 0.6, linewidth = 0.4) +
    geom_hline(yintercept = hist_tmax, linetype = "dashed", colour = "black", linewidth = 0.6) +
    annotate("text", x = -Inf, y = hist_tmax, label = "99% Thermal Max during 1941-2022",
             hjust = -0.05, vjust = -0.4, size = 2.8, colour = "grey20") +
    geom_hline(yintercept = hist_tmin, linetype = "dashed", colour = "black", linewidth = 0.6) +
    annotate("text", x = -Inf, y = hist_tmin, label = "1% Thermal Min during 1941-2022",
             hjust = -0.05, vjust = 1.4, size = 2.8, colour = "grey20") +
    scale_fill_gradientn(
      name    = "Historical CDF (1940-2022)",
      colours = c("#313695", "#74add1", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#d73027"),
      values  = scales::rescale(c(0, cdf_lo, (cdf_lo + 0.5) / 2, 0.5,
                                  (cdf_hi + 0.5) / 2, cdf_hi, 1)),
      limits  = c(0, 1),
      breaks  = c(cdf_lo, 0.25, 0.5, 0.75, cdf_hi),
      labels  = scales::percent_format(accuracy = 1)
    ) +
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top",
                                 barwidth = 15, barheight = 0.5)) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("Distribution of %s across range (%d cells)",
                             VAR_LABELS["temp__12_up"], range_size),
             hjust = -0.03, vjust = 1.5, size = 2.8, colour = "grey30") +
    scale_y_continuous(name = "Annual max temperature (\u00b0C)", expand = c(0.2, 0.2)) +
    scale_x_continuous(name = "Year", breaks = seq(1950, 2025, by = 10), expand = c(0, 0)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )

  # --- Density side panel ----------------------------------------------------
  DENS_ADJUST  <- 0.5
  hist_vals    <- temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY]
  yr_last      <- max(temp_long$year)
  last_vals    <- temp_long$temp_c[temp_long$year == yr_last]

  # Need ≥ 2 distinct points to estimate KDE bandwidth automatically
  if (length(hist_vals) < 2 || length(last_vals) < 2) return(NULL)

  # No from/to clamping — kernel tapers naturally; normalise each curve to max=1
  dens_hist_fn <- density(hist_vals, adjust = DENS_ADJUST, n = N_GRID)
  dens_2025_fn <- density(last_vals, adjust = DENS_ADJUST, n = N_GRID)

  dens_side <- dplyr::bind_rows(
    tibble::tibble(temp_c  = dens_hist_fn$x,
                   density = dens_hist_fn$y / max(dens_hist_fn$y),
                   period  = sprintf("1941\u2013%d", PERIOD_BOUNDARY)),
    tibble::tibble(temp_c  = dens_2025_fn$x,
                   density = dens_2025_fn$y / max(dens_2025_fn$y),
                   period  = as.character(yr_last))
  ) |>
    dplyr::mutate(period = factor(period, levels = c(
      sprintf("1941\u2013%d", PERIOD_BOUNDARY), as.character(yr_last))))

  period_cols <- c("grey40", "#7b2d8b")
  names(period_cols) <- levels(dens_side$period)

  # Align y-axes: density kernel extends slightly beyond data; use its range for
  # both panels so ticks line up perfectly.
  shared_temp_lim <- range(dens_side$temp_c)
  p_kde <- p_kde + ggplot2::scale_y_continuous(
    name   = "Annual max temperature (\u00b0C)",
    limits = shared_temp_lim,
    expand = c(0, 0)
  )

  last_tmax  <- quantile(last_vals, 0.99, type = 7)
  last_tmin  <- quantile(last_vals, 0.01, type = 7)
  delta_hi   <- last_tmax - hist_tmax
  delta_lo   <- last_tmin - hist_tmin
  dens_at_hi <- approx(dens_2025_fn$x, dens_2025_fn$y / max(dens_2025_fn$y),
                       xout = last_tmax)$y
  dens_at_lo <- approx(dens_2025_fn$x, dens_2025_fn$y / max(dens_2025_fn$y),
                       xout = last_tmin)$y

  p_dens <- ggplot(dens_side, aes(x = temp_c, y = density, fill = period, colour = period)) +
    geom_area(alpha = 0.35, position = "identity", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = hist_tmax, linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_vline(xintercept = hist_tmin, linetype = "dashed", colour = "black", linewidth = 0.5) +
    # Arrows + labels showing shift of 99th/1st percentile from hist -> most-recent year
    annotate("segment",
             x = hist_tmax, xend = last_tmax,
             y = dens_at_hi, yend = dens_at_hi,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
             colour = "#7b2d8b", linewidth = 0.7) +
    annotate("text",
             x = (hist_tmax + last_tmax) / 2, y = dens_at_hi + 0.1,
             label = sprintf("%+.1f\u00b0C", delta_hi),
             hjust = 0.5, size = 2.8, colour = "#7b2d8b") +
    annotate("point",
             x = last_tmax, y = dens_at_hi,
             colour = "#7b2d8b", size = 2.5) +
    annotate("segment",
             x = hist_tmin, xend = last_tmin,
             y = dens_at_lo, yend = dens_at_lo,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
             colour = "#7b2d8b", linewidth = 0.7) +
    annotate("text",
             x = (hist_tmin + last_tmin) / 2, y = dens_at_lo + 0.1,
             label = sprintf("%+.1f\u00b0C", delta_lo),
             hjust = 0.5, size = 2.8, colour = "#7b2d8b") +
    annotate("point",
             x = last_tmin, y = dens_at_lo,
             colour = "#7b2d8b", size = 2.5) +
    coord_flip() +
    scale_x_continuous(limits  = shared_temp_lim,
                       expand  = c(0, 0), name = NULL, labels = NULL) +
    scale_y_continuous(name = "Density", expand = c(0, 0)) +
    scale_fill_manual(values = period_cols, name = NULL) +
    scale_colour_manual(values = period_cols, name = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid           = element_blank(),
      axis.ticks.y         = element_blank(),
      legend.position      = c(1, 0.2),
      legend.justification = c(1, 0),
      legend.text          = element_text(size = 8),
      plot.margin          = margin(0, 4, 0, 0)
    )

  # --- Exposure strip (all 8 variables, stacked bars) ------------------------
  all_vars_stacked <- c(
    "temp__12_up", "temp__3__max_up",
    "temp__12_lo", "temp__3__min_lo",
    "precip__12_up", "precip__3__max_up",
    "precip__12_lo", "precip__3__min_lo"
  )
  var_labels_stacked <- c(
    temp__12_up       = "High temperature (annual)",
    temp__3__max_up   = "High temperature (warmest 3 mo)",
    temp__12_lo       = "Low temperature (annual)",
    temp__3__min_lo   = "Low temperature (coldest 3 mo)",
    precip__12_up     = "High precipitation (annual)",
    precip__3__max_up = "High precipitation (wettest 3 mo)",
    precip__12_lo     = "Low precipitation (annual)",
    precip__3__min_lo = "Low precipitation (driest 3 mo)"
  )
  var_colors_stacked <- c(
    "High temperature (annual)"          = "#d73027",
    "High temperature (warmest 3 mo)"    = "#f46d43",
    "Low temperature (annual)"           = "#313695",
    "Low temperature (coldest 3 mo)"     = "#74add1",
    "High precipitation (annual)"        = "#1a9641",
    "High precipitation (wettest 3 mo)"  = "#a6d96a",
    "Low precipitation (annual)"         = "#d9a427",
    "Low precipitation (driest 3 mo)"    = "#fdcc8a"
  )
  recent_year_range <- range(allcell_sp$year)

  accumulation <- compute_exposure_accumulation(allcell_sp, recent_year_range)

  var_exp_stacked <- allcell_sp |>
    dplyr::filter(var %in% all_vars_stacked) |>
    dplyr::group_by(var, year) |>
    dplyr::summarise(n_cells = dplyr::n_distinct(cell), .groups = "drop") |>
    tidyr::complete(
      var  = all_vars_stacked,
      year = seq(recent_year_range[1], recent_year_range[2]),
      fill = list(n_cells = 0L)
    ) |>
    dplyr::mutate(
      pct_exposed = n_cells / range_size * 100,
      var_label   = factor(var_labels_stacked[var],
                           levels = var_labels_stacked[all_vars_stacked])
    )

  p_vars <- ggplot() +
    # Cumulative unique cells ever exposed (any variable) — grey filled area
    geom_area(data = accumulation, aes(x = year, y = cumulative_pct),
              fill = "grey70", colour = NA) +
    # Annual per-variable exposure — stacked bars
    geom_col(data = var_exp_stacked, aes(x = year, y = pct_exposed, fill = var_label),
             width = 0.8) +
    scale_fill_manual(
      values = var_colors_stacked,
      name   = NULL,
      guide  = guide_legend(ncol = 2, byrow = FALSE)
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(1950, 2025, by = 10)) +
    scale_y_continuous(
      name   = "% range exposed",
      labels = function(x) paste0(x, "%"),
      limits = c(0, 100),
      expand = c(0, 0)
    ) +
    labs(x = NULL,
         caption = "Grey area = cumulative unique cells ever exposed; bars = annual per-variable exposure") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      plot.caption       = element_text(size = 7, colour = "grey50"),
      legend.key.size    = unit(0.35, "cm"),
      legend.text        = element_text(size = 7),
      legend.position    = c(0.01, 0.99),
      legend.justification = c(0, 1),
      legend.background  = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.margin      = margin(3, 5, 3, 5)
    )

  # --- Combine ---------------------------------------------------------------
  patchwork::wrap_plots(
    p_vars, patchwork::plot_spacer(), p_kde, p_dens,
    ncol = 2, widths = c(5, 2), heights = c(1.2, 3)
  ) +
    patchwork::plot_layout(axes = "collect_x", guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Historical exposure of %s", gsub("_", " ", sp_name)),
      theme = theme(plot.title    = element_text(face = "italic"),
                    legend.position  = "bottom",
                    legend.direction = "horizontal")
    )
}


# =============================================================================
# 4.  Render all species — write individual PNG files
# =============================================================================

species_figs_dir <- file.path(cache_dir, "species_figs")
dir.create(species_figs_dir, showWarnings = FALSE, recursive = TRUE)

n_rendered <- 0L
n_failed   <- 0L
for (sp_name in all_species) {
  out_path <- file.path(species_figs_dir, paste0(sp_name, "_kde.png"))
  if (!file.exists(out_path)) {
    tryCatch({
      fig <- build_kde_figure(
        sp_name    = sp_name,
        allcell_sp = allcell_list[[sp_name]],
        ts_data    = ts_cache[[sp_name]]
      )
      if (!is.null(fig)) {
        ggplot2::ggsave(out_path, plot = fig, width = 3, height = 2, dpi = 240)
        n_rendered <- n_rendered + 1L
      }
    }, error = function(e) {
      message(sprintf("[kde-cache] Failed for %s: %s", sp_name, e$message))
      n_failed <<- n_failed + 1L
    })
  }
}

cat(sprintf(
  "[kde-cache] Rendered %d figures (%d failed)\n",
  n_rendered, n_failed
))

cat(sprintf("[kde-cache] PNG files written to: %s/\n", species_figs_dir))
cat("[kde-cache] Done.\n")
