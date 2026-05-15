# =============================================================================
# explorefiguresv3.r
#
# Prototype: generalized multi-raster KDE timeseries figure.
# Extends the single-variable KDE (currently only temp__12_up) to show all
# rasters within a variable family in a single faceted figure.
#
# Design: one facet row per raster. Within each row the same KDE-heatmap
# idiom as the current build_kde_figure() applies:
#   - Background tiles  : historical ECDF percentile (blue → yellow → red)
#   - Cell trajectories : grey (normal), colored (exposed)
#   - Threshold line(s) : "up" rasters get an upper dashed line (value > p99)
#                         "lo" rasters get a lower dashed line (value < p01)
#                         rasters with both get both lines (e.g. annual max)
#
# Usage: source() this file interactively with SPECIES set below.
# Output: saves _test_temperature_kde.png
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

source("_R/timeseries-utils.R")
source("_R/timeseries-plotting.R")


# =============================================================================
# 1. Raster configuration
# =============================================================================
# One row per raster file. A single raster can serve both an upper-extreme
# variable (up_var) and a lower-extreme variable (lo_var), or just one.
#
#   "up"  →  threshold = 99th pctile of per-cell 99th pctile (historical)
#            exposure  = value > threshold  (anomalously hot / wet)
#   "lo"  →  threshold = 1st  pctile of per-cell 1st  pctile (historical)
#            exposure  = value < threshold  (anomalously cold / dry)
#
# offset: subtracted from raw raster values (273.15 for Kelvin→°C, 0 for mm).
# Rasters are in data/rast/; adjust RAST_DIR if needed.
# =============================================================================

# NOTE: temp__3__min.tif is currently a byte-for-byte duplicate of
# temp__3__max.tif — the coldest-quarter raster is not yet available.
# Once rebuilt, "Seasonal min" will show the coldest-quarter distribution.
RASTER_CONFIG <- tibble::tribble(
  ~family,         ~raster_file,         ~label,           ~unit,       ~offset, ~up_var,            ~lo_var,            ~color_up,   ~color_lo,
  "Temperature",   "temp__12.tif",       "Annual mean",    "\u00b0C",   273.15,  "temp__12_up",      "temp__12_lo",      "#d73027",   "#313695",
  "Temperature",   "temp__3__max.tif",   "Warmest quarter","\u00b0C",   273.15,  "temp__3__max_up",  NA_character_,      "#f46d43",   NA_character_,
  "Temperature",   "temp__3__min.tif",   "Coldest quarter","\u00b0C",   273.15,  NA_character_,      "temp__3__min_lo",   NA_character_, "#74add1",
  "Precipitation", "precip__12.tif",     "Annual",        "mm",        0,       "precip__12_up",    "precip__12_lo",    "#1a9641",   "#d9a427",
  "Precipitation", "precip__3__max.tif", "Seasonal max",  "mm",        0,       "precip__3__max_up", NA_character_,     "#a6d96a",   NA_character_,
  "Precipitation", "precip__3__min.tif", "Seasonal min",  "mm",        0,       NA_character_,      "precip__3__min_lo", NA_character_, "#fdcc8a"
)


# =============================================================================
# 2. Test configuration — adjust to explore different species
# =============================================================================

DATA_DIR        <- "data"
RAST_DIR        <- file.path(DATA_DIR, "rast")
ALLCELL         <- file.path(DATA_DIR, "AllCellExposureSpXVar.rds")
SPECIES         <- "Hynobius_retardatus"
N_GRID          <- 300L
PERIOD_BOUNDARY <- 2022L
N_CELL_LINES    <- 200L    # subsampled cells for trajectory lines


# =============================================================================
# 3. Load species data
# =============================================================================

cat("Loading allcell data...\n")
allcell_all <- readRDS(ALLCELL)
allcell_sp  <- dplyr::filter(allcell_all, spName == SPECIES)
range_cells <- unique(allcell_sp$cell)
range_size  <- allcell_sp$rangeSize[1]
cat(sprintf("Species: %s | range: %d cells\n", SPECIES, range_size))
rm(allcell_all); gc()


# =============================================================================
# 4. Helper: extract_raster_values()
#
# Loads one raster, extracts per-cell time series for the species' range cells,
# converts units, and returns a long-format tibble (cell, year, value).
# Returns NULL with a message if the raster file is not found.
# =============================================================================

extract_raster_values <- function(raster_file, range_cells, unit_offset,
                                  rast_dir = RAST_DIR) {
  r_path <- file.path(rast_dir, raster_file)
  if (!file.exists(r_path)) {
    message(sprintf("  [skip] Raster not found: %s", r_path))
    return(NULL)
  }
  r           <- terra::rast(r_path)
  layer_years <- as.integer(sub("WY", "", names(r)))
  raw_mat     <- terra::extract(r, range_cells)

  as.data.frame(raw_mat) |>
    dplyr::mutate(cell = range_cells) |>
    tidyr::pivot_longer(cols = -cell, names_to = "layer", values_to = "raw_val") |>
    dplyr::mutate(
      year  = layer_years[match(layer, names(r))],
      value = raw_val - unit_offset
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::select(cell, year, value)
}


# =============================================================================
# 5. Helper: compute_var_kde_data()
#
# Builds the KDE tile grid for one raster's values:
#   - Fits a historical ECDF from all cell-year values up to period_boundary
#   - Evaluates that ECDF on a fine grid (n_grid points) for every year
#   - Clips tiles to the per-year observed envelope (no overhang)
#
# Returns: list(tiles, hist_ecdf, range_lines)
# =============================================================================

compute_var_kde_data <- function(val_long, n_grid = N_GRID,
                                 period_boundary = PERIOD_BOUNDARY) {
  range_lines <- val_long |>
    dplyr::group_by(year) |>
    dplyr::summarise(val_min = min(value), val_max = max(value), .groups = "drop")

  val_range <- range(val_long$value)
  grid_vals <- seq(val_range[1], val_range[2], length.out = n_grid)
  grid_step <- diff(grid_vals)[1]
  all_years <- sort(unique(val_long$year))

  hist_ecdf <- stats::ecdf(val_long$value[val_long$year <= period_boundary])

  tiles <- tidyr::expand_grid(year = all_years, val_mid = grid_vals) |>
    dplyr::mutate(pct_range = hist_ecdf(val_mid)) |>
    dplyr::left_join(range_lines, by = "year") |>
    dplyr::filter(val_mid >= val_min - grid_step / 2,
                  val_mid <= val_max + grid_step / 2) |>
    dplyr::mutate(
      ymin = pmax(val_mid - grid_step / 2, val_min),
      ymax = pmin(val_mid + grid_step / 2, val_max)
    )

  list(tiles = tiles, hist_ecdf = hist_ecdf, range_lines = range_lines)
}


# =============================================================================
# 6. Helper: compute_threshold()
#
# Computes the exposure threshold for one raster × direction:
#   "up" → 99th pctile of per-cell 99th pctile  (upper extreme)
#   "lo" →  1st pctile of per-cell 1st  pctile  (lower extreme)
#
# Only historical values (year <= period_boundary) are used.
# =============================================================================

compute_threshold <- function(val_long, direction,
                              period_boundary = PERIOD_BOUNDARY) {
  hist <- dplyr::filter(val_long, year <= period_boundary)
  q    <- if (direction == "up") 0.99 else 0.01

  hist |>
    dplyr::group_by(cell) |>
    dplyr::summarise(cell_extreme = quantile(value, q), .groups = "drop") |>
    dplyr::summarise(val = quantile(cell_extreme, q)) |>
    dplyr::pull(val)
}


# =============================================================================
# 7. Main: build_family_kde_figure()
#
# For a given variable family (e.g. "Temperature"), loops over all rasters in
# raster_config and assembles a two-panel patchwork:
#
#   Left  (wide)  — KDE heatmap: one facet row per raster, free y-scale,
#                   background tiles colored by historical CDF percentile,
#                   individual cell trajectories, and threshold dashed line(s).
#   Right (narrow) — Density side panel: historical vs. most-recent-year KDE
#                    curves per facet row, with matching threshold vlines.
#
# Arguments:
#   family_name   – character; must match a value in raster_config$family
#   raster_config – RASTER_CONFIG tibble defined above
#   range_cells   – integer vector of terra cell indices for the focal species
#   n_grid        – KDE grid resolution (default N_GRID)
#   period_boundary – last year of the historical reference period
#   n_cell_lines  – number of cells to subsample for trajectory lines
#   rast_dir      – directory containing the raster .tif files
# =============================================================================

build_family_kde_figure <- function(family_name, raster_config, range_cells,
                                    n_grid          = N_GRID,
                                    period_boundary = PERIOD_BOUNDARY,
                                    n_cell_lines    = N_CELL_LINES,
                                    rast_dir        = RAST_DIR) {

  fam_cfg    <- dplyr::filter(raster_config, family == family_name)
  unit_label <- fam_cfg$unit[1]
  cat(sprintf("[%s] Processing %d rasters...\n", family_name, nrow(fam_cfg)))


  # --- 1. Extract raster values for each config row ------------------------
  rast_data <- list()
  for (i in seq_len(nrow(fam_cfg))) {
    cfg  <- fam_cfg[i, ]
    vals <- extract_raster_values(cfg$raster_file, range_cells, cfg$offset, rast_dir)
    if (!is.null(vals)) rast_data[[cfg$raster_file]] <- vals
  }

  present <- names(rast_data)
  if (length(present) == 0) {
    message(sprintf("[%s] No rasters found — skipping.", family_name))
    return(NULL)
  }
  cat(sprintf("[%s] Loaded %d/%d rasters.\n", family_name, length(present), nrow(fam_cfg)))


  # --- 2. KDE tiles and thresholds per raster ------------------------------
  kde_list    <- list()
  range_list  <- list()   # per-year min/max envelope across all range cells
  thresh_rows <- list()

  for (rfile in present) {
    cfg     <- dplyr::filter(fam_cfg, raster_file == rfile)
    vals    <- rast_data[[rfile]]
    kde_res <- compute_var_kde_data(vals, n_grid, period_boundary)

    kde_list[[rfile]] <- kde_res$tiles |>
      dplyr::mutate(raster_id = cfg$label)

    range_list[[rfile]] <- kde_res$range_lines |>
      dplyr::mutate(raster_id = cfg$label)

    # Upper threshold row (if raster has an up_var)
    if (!is.na(cfg$up_var)) {
      t_up <- compute_threshold(vals, "up", period_boundary)
      thresh_rows[[paste0(rfile, "_up")]] <- tibble::tibble(
        raster_id  = cfg$label,
        threshold  = t_up,
        direction  = "up",
        color      = cfg$color_up,
        thresh_lbl = sprintf("99th pctile upper (%.1f %s)", t_up, unit_label),
        text_vjust = -0.4
      )
      cat(sprintf("  %s [up]: threshold = %.2f %s\n", cfg$label, t_up, unit_label))
    }

    # Lower threshold row (if raster has a lo_var)
    if (!is.na(cfg$lo_var)) {
      t_lo <- compute_threshold(vals, "lo", period_boundary)
      thresh_rows[[paste0(rfile, "_lo")]] <- tibble::tibble(
        raster_id  = cfg$label,
        threshold  = t_lo,
        direction  = "lo",
        color      = cfg$color_lo,
        thresh_lbl = sprintf("1st pctile lower (%.1f %s)", t_lo, unit_label),
        text_vjust = 1.4
      )
      cat(sprintf("  %s [lo]: threshold = %.2f %s\n", cfg$label, t_lo, unit_label))
    }
  }


  # --- 3. Subsampled cell trajectory segments per raster -------------------
  set.seed(42)
  cell_sample <- sample(range_cells, min(n_cell_lines, length(range_cells)))
  seg_list    <- list()

  for (rfile in present) {
    cfg   <- dplyr::filter(fam_cfg, raster_file == rfile)
    vals  <- dplyr::filter(rast_data[[rfile]], cell %in% cell_sample)
    t_up  <- if (!is.na(cfg$up_var)) thresh_rows[[paste0(rfile, "_up")]]$threshold else NA_real_
    t_lo  <- if (!is.na(cfg$lo_var)) thresh_rows[[paste0(rfile, "_lo")]]$threshold else NA_real_

    segs <- vals |>
      dplyr::arrange(cell, year) |>
      dplyr::group_by(cell) |>
      dplyr::mutate(
        year_end   = dplyr::lead(year),
        val_end    = dplyr::lead(value),
        exposed_up = !is.na(t_up) & value > t_up,
        exposed_lo = !is.na(t_lo) & value < t_lo,
        # Pre-compute display color; upper takes priority if both triggered
        exp_color  = dplyr::case_when(
          exposed_up ~ cfg$color_up,
          exposed_lo ~ cfg$color_lo,
          TRUE       ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(year_end)) |>
      dplyr::ungroup() |>
      dplyr::mutate(raster_id = cfg$label)

    seg_list[[rfile]] <- segs
  }


  # --- 4. Combine into faceted data frames ---------------------------------
  label_order  <- dplyr::filter(fam_cfg, raster_file %in% present)$label

  tiles_df    <- dplyr::bind_rows(kde_list) |>
    dplyr::mutate(raster_id = factor(raster_id, levels = label_order))
  envelope_df <- dplyr::bind_rows(range_list) |>
    dplyr::mutate(raster_id = factor(raster_id, levels = label_order))
  segs_df     <- dplyr::bind_rows(seg_list) |>
    dplyr::mutate(raster_id = factor(raster_id, levels = label_order))
  thresh_df   <- dplyr::bind_rows(thresh_rows) |>
    dplyr::mutate(raster_id = factor(raster_id, levels = label_order))

  segs_normal  <- dplyr::filter(segs_df,  is.na(exp_color))
  segs_exposed <- dplyr::filter(segs_df, !is.na(exp_color))

  cdf_cols <- c("#313695", "#4575b4", "#74add1", "#e0f3f8",
                "#ffffbf", "#fee090", "#fdae61", "#d73027")


  # --- 5. KDE heatmap panel ------------------------------------------------
  p_kde <- ggplot(tiles_df, aes(fill = pct_range)) +
    geom_rect(aes(xmin = year - 0.5, xmax = year + 0.5,
                  ymin = ymin, ymax = ymax)) +

    # Unexposed cell trajectories — thin grey lines
    geom_segment(
      data        = segs_normal,
      aes(x = year, xend = year_end, y = value, yend = val_end),
      inherit.aes = FALSE,
      colour = "black", alpha = 0.10, linewidth = 0.20
    ) +

    # Exposed cell trajectories — colored by pre-computed exposure color
    geom_segment(
      data        = segs_exposed,
      aes(x = year, xend = year_end, y = value, yend = val_end,
          colour = I(exp_color)),
      inherit.aes = FALSE,
      alpha = 0.55, linewidth = 0.35
    ) +

    # Range envelope: thick step lines bounding the min/max across all cells.
    # Shift x by -0.5 so each horizontal segment spans year-0.5 to year+0.5,
    # matching the geom_rect tiles which are centered on each year.
    geom_step(
      data        = envelope_df,
      aes(x = year - 0.5, y = val_max),
      inherit.aes = FALSE,
      colour = "black", linewidth = 0.8
    ) +
    geom_step(
      data        = envelope_df,
      aes(x = year - 0.5, y = val_min),
      inherit.aes = FALSE,
      colour = "black", linewidth = 0.8
    ) +

    # Threshold dashed line(s) per facet — color matches variable direction
    geom_hline(
      data        = thresh_df,
      aes(yintercept = threshold, colour = I(color)),
      linetype    = "dashed", linewidth = 0.55,
      inherit.aes = FALSE
    ) +

    # Threshold label: "up" labels sit above the line, "lo" labels below
    geom_text(
      data        = thresh_df,
      aes(x = -Inf, y = threshold, label = thresh_lbl,
          vjust = text_vjust),
      hjust       = -0.03, size = 2.6, colour = "grey25",
      inherit.aes = FALSE
    ) +

    facet_grid(rows = vars(raster_id), scales = "free_y") +

    scale_fill_gradientn(
      name    = sprintf("Historical CDF (1941\u2013%d)", period_boundary),
      colours = cdf_cols,
      limits  = c(0, 1),
      labels  = scales::percent_format(accuracy = 1)
    ) +
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top",
                                 barwidth = 13, barheight = 0.5)) +
    scale_x_continuous(breaks = seq(1950, 2030, by = 10), expand = c(0, 0)) +
    scale_y_continuous(
      name   = paste0("Value (", unit_label, ")"),
      expand = c(0.12, 0.12)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.background   = element_rect(fill = "grey93", colour = "grey72"),
      strip.text         = element_text(size = 9, face = "bold"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )


  # --- 6. Density side panel -----------------------------------------------
  # Use the minimum of the per-raster last years so every raster has data for
  # the "recent" comparison year (annual raster may run one year ahead of the
  # seasonal rasters).
  DENS_ADJUST  <- 0.5
  dens_rows    <- list()
  last_yr_lbl  <- min(vapply(rast_data, function(d) max(d$year), integer(1)))

  for (rfile in present) {
    cfg          <- dplyr::filter(fam_cfg, raster_file == rfile)
    vals         <- rast_data[[rfile]]
    hist_vals    <- vals$value[vals$year <= period_boundary]
    # Use this raster's last available year (not the global max) to avoid
    # empty recent_vals when rasters cover different year ranges.
    rast_last_yr <- max(vals$year)
    recent_vals  <- vals$value[vals$year == rast_last_yr]
    if (length(hist_vals) < 2 || length(recent_vals) < 2) next

    dh <- density(hist_vals,   adjust = DENS_ADJUST, n = 256)
    dr <- density(recent_vals, adjust = DENS_ADJUST, n = 256)

    dens_rows[[rfile]] <- dplyr::bind_rows(
      tibble::tibble(
        value     = dh$x,
        density   = dh$y / max(dh$y),
        period    = sprintf("1941\u2013%d", period_boundary),
        raster_id = cfg$label
      ),
      tibble::tibble(
        value     = dr$x,
        density   = dr$y / max(dr$y),
        period    = as.character(last_yr_lbl),  # consistent display label
        raster_id = cfg$label
      )
    )
  }

  # Add invisible anchor points at the per-row KDE y-range extremes so the
  # density panel x-scale (which becomes y after coord_flip) aligns with the
  # KDE heatmap y-scale for each facet row.
  row_y_ranges <- tiles_df |>
    dplyr::group_by(raster_id) |>
    dplyr::summarise(v_lo = min(ymin), v_hi = max(ymax), .groups = "drop")

  anchor_rows <- dplyr::bind_rows(
    dplyr::transmute(row_y_ranges, raster_id = as.character(raster_id),
                     value = v_lo, density = 0, period = NA_character_),
    dplyr::transmute(row_y_ranges, raster_id = as.character(raster_id),
                     value = v_hi, density = 0, period = NA_character_)
  )

  dens_df <- dplyr::bind_rows(dens_rows) |>
    dplyr::bind_rows(anchor_rows) |>
    dplyr::mutate(
      raster_id = factor(raster_id, levels = label_order),
      period    = factor(period, levels = c(
        sprintf("1941\u2013%d", period_boundary), as.character(last_yr_lbl)
      ))
    )

  period_cols <- setNames(
    c("grey45", "#7b2d8b"),
    c(sprintf("1941\u2013%d", period_boundary), as.character(last_yr_lbl))
  )

  p_dens <- ggplot(dens_df, aes(x = value, y = density,
                                colour = period, fill = period)) +
    geom_area(alpha = 0.32, position = "identity", linewidth = 0.55) +

    # Threshold vlines matching the heatmap dashes (before coord_flip,
    # these vertical lines become horizontal after flipping — matching the
    # dashed hlines in p_kde)
    geom_vline(
      data        = thresh_df,
      aes(xintercept = threshold, colour = I(color)),
      linetype    = "dashed", linewidth = 0.45,
      inherit.aes = FALSE
    ) +

    facet_grid(rows = vars(raster_id), scales = "free") +
    coord_flip() +
    scale_colour_manual(values = period_cols, name = NULL, na.translate = FALSE) +
    scale_fill_manual(values   = period_cols, name = NULL, na.translate = FALSE) +
    scale_x_continuous(labels = NULL, name = NULL, expand = c(0, 0)) +
    scale_y_continuous(name = "Density", expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text      = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )


  # --- 7. Combine KDE + density panels -------------------------------------
  p_kde + p_dens +
    patchwork::plot_layout(widths = c(4, 1), guides = "collect") +
    patchwork::plot_annotation(
      title    = sprintf("%s \u2014 %s distribution KDE",
                         gsub("_", " ", SPECIES), family_name),
      subtitle = sprintf(
        "Background = historical CDF percentile (1941\u2013%d)  |  Range: %d cells  |  %d cell trajectories subsampled",
        period_boundary, length(range_cells), min(n_cell_lines, length(range_cells))
      ),
      theme = theme(
        plot.title       = element_text(face = "italic", size = 13),
        plot.subtitle    = element_text(size = 8.5, colour = "grey45"),
        legend.position  = "bottom",
        legend.direction = "horizontal"
      )
    )
}


# =============================================================================
# 8. Test: Temperature family — single species, single figure
# =============================================================================

cat("\n--- Test: Temperature family ---\n")
fig_temp <- build_family_kde_figure(
  family_name   = "Temperature",
  raster_config = RASTER_CONFIG,
  range_cells   = range_cells
)

if (!is.null(fig_temp)) {
  out_path <- "_test_temperature_kde.png"
  ggplot2::ggsave(out_path, fig_temp,
                  width = 14, height = 11, dpi = 110, bg = "white")
  cat(sprintf("Saved: %s\n", out_path))
}


# =============================================================================
# Precipitation test — uncomment once rasters are confirmed
# =============================================================================
# cat("\n--- Test: Precipitation family ---\n")
# fig_precip <- build_family_kde_figure(
#   family_name   = "Precipitation",
#   raster_config = RASTER_CONFIG,
#   range_cells   = range_cells
# )
# if (!is.null(fig_precip)) {
#   ggplot2::ggsave("_test_precipitation_kde.png", fig_precip,
#                   width = 14, height = 11, dpi = 110, bg = "white")
#   cat("Saved: _test_precipitation_kde.png\n")
# }
