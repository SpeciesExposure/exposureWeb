# =============================================================================
# 02_explore_figures_v2.R
#
# Exploratory composite figure for a single species showing:
#   (a) p_vars    — stacked bar: annual % of range exposed per variable
#   (b) p_map     — species range map coloured by 2025 CDF position,
#                   with CartoDB basemap and orthographic globe inset
#   (c) p         — temperature heatmap (year × temp, coloured by hist CDF)
#   (d) p_dens    — density curves (historical vs most-recent year) with
#                   percentile shift annotations
#
# The four panels are assembled into p_combined (2×2 patchwork) and also
# produced as three standalone composites: p_combined_vars, p_combined_ts,
# p_combined_map.  A polar summary chart (p_polar) is produced at the end.
#
# Data sources:
#   data/temp__12.tif              – annual max temp, WY1941–WY2025
#   data/AllCellExposureSpXVar.rds – per-cell exposure flags by variable
# =============================================================================

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

source("_R/timeseries-utils.R")     # classify_time_period

# --- Config -------------------------------------------------------------------

SPECIES    <- "Hynobius_retardatus"  # change to explore any species
RASTER     <- "data/temp__12.tif"    # annual max temperature raster (Kelvin)
ALLCELL    <- "data/AllCellExposureSpXVar.rds"  # exposure flags by cell × var × year
N_GRID     <- 1000L                  # KDE / CDF evaluation grid resolution

# --- Load species exposure records -------------------------------------------
# One row per (cell, var, year) where the cell exceeded the exposure threshold.

cat("Reading allcell data...\n")
allcell <- readRDS(ALLCELL) |>
  dplyr::filter(spName == SPECIES)

range_cells <- unique(allcell$cell)   # terra cell indices into the 588×1440 grid
range_size  <- allcell$rangeSize[1]   # total occupied cells (denominator for %)
cat(sprintf("Species: %s | range size: %d | observed cells: %d\n",
            SPECIES, range_size, length(range_cells)))

# --- Extract temperature time series for range cells -------------------------

cat("Reading temperature raster...\n")
r <- terra::rast(RASTER)

# Layer names are "WY1941" … "WY2025"; strip prefix to get integer years
layer_years <- as.integer(sub("WY", "", names(r)))

# Extract all layers at once for the range cells (rows = cells, cols = years)
cat("Extracting values for range cells...\n")
temp_mat <- terra::extract(r, range_cells)

# Convert Kelvin -> Celsius and reshape to long format
temp_long <- temp_mat |>
  as.data.frame() |>
  dplyr::mutate(cell = range_cells) |>
  tidyr::pivot_longer(
    cols      = -cell,
    names_to  = "layer",
    values_to = "temp_k"
  ) |>
  dplyr::mutate(
    year   = layer_years[match(layer, names(r))],
    temp_c = temp_k - 273.15
  ) |>
  dplyr::filter(!is.na(temp_c), year >= 1941, year <= 2025)

cat(sprintf("Temperature range: %.1f – %.1f °C\n",
            min(temp_long$temp_c), max(temp_long$temp_c)))

# --- Per-year temperature envelope and CDF-coloured grid --------------------
# We evaluate the historical CDF on a shared fine temperature grid so every
# year's geom_rect tiles have identical height, producing a smooth gradient.

# Per-year observed min/max — used to clip tiles to the actual data envelope
temp_range_lines <- temp_long |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    temp_min = min(temp_c),
    temp_max = max(temp_c),
    .groups  = "drop"
  )

temp_range  <- range(temp_long$temp_c)
grid_temps  <- seq(temp_range[1], temp_range[2], length.out = N_GRID)
grid_step   <- diff(grid_temps)[1]

all_years_v <- sort(unique(temp_long$year))

# Historical reference period: CDF is built from 1941–PERIOD_BOUNDARY only
PERIOD_BOUNDARY <- 2022L

# Historical CDF built once from the reference period; applied to all years
# so the fill colour always represents position within the historical baseline.
hist_ecdf_fn  <- stats::ecdf(temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY])

# Tile grid: for each (year, temp grid point) inside the observed envelope,
# record the historical CDF value. ymin/ymax are clipped to the actual
# per-year min/max so tiles don't overhang the envelope lines.
temp_binned <- tidyr::expand_grid(
  year     = all_years_v,
  temp_mid = grid_temps
) |>
  dplyr::mutate(pct_range = hist_ecdf_fn(temp_mid)) |>
  dplyr::left_join(temp_range_lines, by = "year") |>
  dplyr::filter(temp_mid >= temp_min - grid_step / 2,
                temp_mid <= temp_max + grid_step / 2) |>
  dplyr::mutate(
    ymin = pmax(temp_mid - grid_step / 2, temp_min),
    ymax = pmin(temp_mid + grid_step / 2, temp_max)
  )

# Exposure thresholds: 99th / 1st percentile of per-cell extremes during
# the historical period.  A cell-year is "exposed" when it exceeds these.
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

# CDF values at the exposure thresholds — used as gradient breakpoints
cdf_lo <- hist_ecdf_fn(hist_tmin)
cdf_hi <- hist_ecdf_fn(hist_tmax)

# Historical median temperature — used for the dashed median line on the heatmap
hist_tmed <- median(temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY])

# --- Heatmap: temperature distribution across range cells over time -----------

# Subsample cells for individual trajectory lines to avoid overplotting
N_CELL_LINES <- 300L
set.seed(42)
cell_sample <- sample(range_cells, min(N_CELL_LINES, length(range_cells)))

# Per-segment data: each row is one year-to-next-year step for a cell.
# Segments are coloured red (exposed hot), blue (exposed cold), or grey.
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

p <- ggplot(temp_binned, aes(fill = pct_range)) +
  geom_rect(aes(xmin = year - 0.5, xmax = year + 0.5, ymin = ymin, ymax = ymax)) +
  geom_line(data = temp_range_lines, aes(x = year, y = temp_min),
            inherit.aes = FALSE, colour = "black", linewidth = 1) +
  geom_line(data = temp_range_lines, aes(x = year, y = temp_max),
            inherit.aes = FALSE, colour = "black", linewidth = 1) +

  # Within-range cell trajectories grey; exposed segments coloured
  geom_segment(data = dplyr::filter(cell_segments, !exposed_hi, !exposed_lo),
               aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
               inherit.aes = FALSE, colour = "black", alpha = 0.12, linewidth = 0.25) +
  geom_segment(data = dplyr::filter(cell_segments, exposed_hi),
               aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
               inherit.aes = FALSE, colour = "#d73027", alpha = 0.6, linewidth = 0.4) +
  geom_segment(data = dplyr::filter(cell_segments, exposed_lo),
               aes(x = year, xend = year_end, y = temp_c, yend = temp_end),
               inherit.aes = FALSE, colour = "#313695", alpha = 0.6, linewidth = 0.4) +
  geom_hline(yintercept = hist_tmax,
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = -Inf, y = hist_tmax, label = "99% Thermal Max",
           hjust = -0.05, vjust = -0.4, size = 2.8, colour = "grey20") +
  geom_hline(yintercept = hist_tmed,
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  geom_hline(yintercept = hist_tmin,
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = -Inf, y = hist_tmin, label = "1% Thermal Min",
           hjust = -0.05, vjust = 1.4, size = 2.8, colour = "grey20") +
  scale_fill_gradientn(
    name    = "Historical CDF (1940-2022)",
    colours = c("#313695", "#74add1", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#d73027"),
    values  = scales::rescale(c(0, cdf_lo, (cdf_lo + 0.5) / 2, 0.5, (cdf_hi + 0.5) / 2, cdf_hi, 1)),
    limits  = c(0, 1),
    breaks  = c(cdf_lo, 0.25, 0.5, 0.75, cdf_hi),
    labels  = scales::percent_format(accuracy = 1)
  ) +
  guides(fill = guide_colorbar(direction = "horizontal", title.position = "top",
                               barwidth = 15, barheight = 0.5)) +
  scale_y_continuous(name = "Annual max temperature (\u00b0C)", expand = c(0.2, 0.2)) +
  scale_x_continuous(name   = "Year",
                     breaks = seq(1950, 2025, by = 10), expand = c(0, 0)) +
  annotate("text",
           x = -Inf, y = Inf,
           label = sprintf("Distribution of %s across range (%d cells)",
                           "Temperature: annual upper extreme", range_size),
           hjust = -0.03, vjust = 1.5, size = 3.2, colour = "grey30") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x           = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.x    = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.major.y    = element_blank(),
    panel.grid.minor      = element_blank(),
    legend.position       = "bottom"
  )
p

# --- Density panel (right side of heatmap) -----------------------------------
# Overlapping area curves for historical baseline vs most-recent year.
# x-axis = density (coord_flip makes it vertical, aligning with heatmap y).

DENS_ADJUST  <- 0.5   # KDE bandwidth multiplier (< 1 = rougher, > 1 = smoother)

hist_vals <- temp_long$temp_c[temp_long$year <= PERIOD_BOUNDARY]
yr_last   <- max(temp_long$year)
last_vals <- temp_long$temp_c[temp_long$year == yr_last]

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
    sprintf("1941\u2013%d", PERIOD_BOUNDARY), as.character(max(temp_long$year)))))

period_cols <- c("grey40", "#7b2d8b")
names(period_cols) <- levels(dens_side$period)

# Align y-axes: density kernel extends slightly beyond data; use its range for
# both the density panel and the heatmap so ticks line up perfectly.
shared_temp_lim <- range(dens_side$temp_c)
p <- p + ggplot2::scale_y_continuous(
  name   = "Annual max temperature (\u00b0C)",
  limits = shared_temp_lim,
  expand = c(0, 0)
)

p_dens <- ggplot(dens_side, aes(x = temp_c, y = density,
                                 fill = period, colour = period)) +
  geom_area(alpha = 0.35, position = "identity", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  # Historical 99th / 1st percentile dashed lines
  geom_vline(xintercept = hist_tmax,
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_vline(xintercept = hist_tmin,
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  # Arrows + labels showing shift of 99th/1st percentile from hist → most-recent year
  {
    last_tmax  <- quantile(last_vals, 0.99, type = 7)
    last_tmin  <- quantile(last_vals, 0.01, type = 7)
    delta_hi   <- last_tmax - hist_tmax
    delta_lo   <- last_tmin - hist_tmin
    dens_at_hi <- approx(dens_2025_fn$x, dens_2025_fn$y / max(dens_2025_fn$y),
                         xout = last_tmax)$y
    dens_at_lo <- approx(dens_2025_fn$x, dens_2025_fn$y / max(dens_2025_fn$y),
                         xout = last_tmin)$y
    list(
      annotate("segment",
               x = hist_tmax, xend = last_tmax,
               y = dens_at_hi, yend = dens_at_hi,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "#7b2d8b", linewidth = 0.7),
      annotate("text",
               x = (hist_tmax + last_tmax) / 2, y = dens_at_hi+.1,
               label = sprintf("%+.1f\u00b0C", delta_hi),
               hjust = 0.5, size = 2.8, colour = "#7b2d8b"),
      annotate("point",
               x = last_tmax, y = dens_at_hi,
               colour = "#7b2d8b", size = 2.5),
      annotate("segment",
               x = hist_tmin, xend = last_tmin,
               y = dens_at_lo, yend = dens_at_lo,
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               colour = "#7b2d8b", linewidth = 0.7),
      annotate("text",
               x = (hist_tmin + last_tmin) / 2, y = dens_at_lo+.1,
               label = sprintf("%+.1f\u00b0C", delta_lo),
               hjust = 0.5, size = 2.8, colour = "#7b2d8b"),
      annotate("point",
               x = last_tmin, y = dens_at_lo,
               colour = "#7b2d8b", size = 2.5)
    )
  } +
  coord_flip() +
  scale_x_continuous(limits  = shared_temp_lim,
                     expand  = c(0, 0),
                     name    = NULL,
                     labels  = NULL) +
  scale_y_continuous(name = "Density", expand = c(0, 0)) +
  scale_fill_manual(values = period_cols, name = NULL) +
  scale_colour_manual(values = period_cols, name = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid       = element_blank(),
    axis.ticks.y     = element_blank(),
    legend.position  = c(1, 0.2),
    legend.justification = c(1, 0),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(0, 4, 0, 0)
  )

# Full allcell table is used directly; alias + year range for convenience
allcell_recent    <- allcell
recent_year_range <- range(allcell_recent$year)

# Cumulative unique cells ever exposed across all variables and years
accumulation <- compute_exposure_accumulation(allcell_recent, recent_year_range)


# Per-variable annual % of range exposed (cells may overlap across variables).
# Grey filled area = cumulative unique cells ever exposed by any variable.

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

var_exp_stacked <- allcell_recent |>
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
    expand = expansion(mult = c(0, 0.08))
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

# --- Map panel (p_map) -------------------------------------------------------
# Range cells coloured by each cell's most-recent-year temperature expressed
# as its percentile in the historical CDF.  CartoDB basemap underneath.

# Per-cell most-recent-year temperature → historical CDF position
cell_hist_map <- temp_long |>
  dplyr::filter(year == max(temp_long$year)) |>
  dplyr::select(cell, temp_c) |>
  dplyr::mutate(cdf_val = hist_ecdf_fn(temp_c))

# XY coordinates for range cells
range_xy_map <- terra::xyFromCell(r[[1]], range_cells) |>
  as.data.frame() |>
  dplyr::mutate(cell = range_cells) |>
  dplyr::left_join(cell_hist_map, by = "cell")

# Map extent: range bounding box + 3° buffer, then xlim widened so the map
# fills the top-right panel of p_combined (14×9", widths=c(5,2), heights=c(1.2,3)).
buf      <- 3
map_xlim <- range(range_xy_map$x) + c(-buf, buf)
map_ylim <- range(range_xy_map$y) + c(-buf, buf)

# Expand xlim to match the panel's geographic aspect ratio
{
  panel_ar  <- (2/7 * 14) / (1.2/4.2 * 9)          # panel width / height in inches
  ctr_lon   <- mean(map_xlim)
  ctr_lat   <- mean(map_ylim)
  lat_span  <- diff(map_ylim)
  lon_span  <- lat_span * panel_ar / cos(ctr_lat * pi / 180)
  map_xlim  <- ctr_lon + c(-0.5, 0.5) * lon_span
}

# Fetch CartoDB.Positron basemap tiles for the map extent
{
  library(maptiles)
  library(tidyterra)
  map_tile_srs <- terra::rast(
    ext = terra::ext(map_xlim[1], map_xlim[2], map_ylim[1], map_ylim[2]),
    crs = "EPSG:4326", nrows = 1L, ncols = 1L
  )
  map_tiles_raw <- maptiles::get_tiles(map_tile_srs,
                                       provider = "CartoDB.Positron",
                                       zoom = 7, crop = TRUE)
  map_tiles <- terra::project(map_tiles_raw, "EPSG:4326", method = "near")
}

p_map <- ggplot() +
  tidyterra::geom_spatraster_rgb(data = map_tiles, maxcell = Inf) +
  geom_tile(data = range_xy_map, aes(x = x, y = y, fill = cdf_val),
            colour = "black", linewidth = 0.15) +
  scale_fill_gradientn(
    colours = c("#313695", "#74add1", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#d73027"),
    values  = scales::rescale(c(0, cdf_lo, (cdf_lo + 0.5) / 2, 0.5, (cdf_hi + 0.5) / 2, cdf_hi, 1)),
    limits  = c(0, 1),
    guide   = "none"
  ) +
  coord_sf(crs = 4326, xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  theme_void()

# --- Globe inset (p_globe) ---------------------------------------------------
# Orthographic projection centred on the species range; hemisphere pre-filtered
# in WGS84 to avoid topology errors when reprojecting antipodal geometries.
library(sf)

globe_lon <- mean(range(range_xy_map$x))
globe_lat <- mean(range(range_xy_map$y))
ortho_crs <- sprintf("+proj=ortho +lat_0=%.2f +lon_0=%.2f +datum=WGS84",
                     globe_lat, globe_lon)

world_sf   <- sf::st_as_sf(maps::map("world", fill = TRUE, plot = FALSE)) |>
  sf::st_make_valid()

# Keep only polygons whose centroid is within ~89° of the focus (visible hemisphere)
world_ctrs <- suppressWarnings(sf::st_centroid(world_sf))
focus_pt   <- sf::st_sfc(sf::st_point(c(globe_lon, globe_lat)), crs = sf::st_crs(world_ctrs))
dists_deg  <- as.numeric(sf::st_distance(world_ctrs, focus_pt)) / 111320
world_hemi <- world_sf[dists_deg < 89, ] |> sf::st_make_valid()

world_ortho <- suppressWarnings(sf::st_transform(world_hemi, ortho_crs)) |>
  sf::st_make_valid()

earth_r <- 6371000
hemi_sfc <- sf::st_sfc(sf::st_point(c(0, 0)), crs = ortho_crs) |>
  sf::st_buffer(dist = earth_r * 0.9999)

# Densify the range bounding box edges before projecting to avoid great-circle
# distortion that would otherwise make the red rectangle look curved.
densify_edge <- function(x0, y0, x1, y1, n = 50) {
  cbind(seq(x0, x1, length.out = n), seq(y0, y1, length.out = n))
}
bbox_pts <- rbind(
  densify_edge(map_xlim[1], map_ylim[1], map_xlim[2], map_ylim[1]),
  densify_edge(map_xlim[2], map_ylim[1], map_xlim[2], map_ylim[2]),
  densify_edge(map_xlim[2], map_ylim[2], map_xlim[1], map_ylim[2]),
  densify_edge(map_xlim[1], map_ylim[2], map_xlim[1], map_ylim[1])
)
bbox_sf <- sf::st_sfc(sf::st_polygon(list(bbox_pts)), crs = 4326) |>
  sf::st_transform(ortho_crs)

p_globe <- ggplot() +
  geom_sf(data = hemi_sfc,   fill = "#ddeeff", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = world_ortho, fill = "grey70", colour = NA) +
  geom_sf(data = bbox_sf,    fill = NA, colour = "red", linewidth = 1.2) +
  # Clip to hemisphere circle via coord_sf limits (no st_intersection needed)
  coord_sf(xlim = c(-earth_r, earth_r), ylim = c(-earth_r, earth_r),
           crs = ortho_crs, expand = FALSE) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.tag         = element_blank()
  )

# Embed globe as inset in the upper-right corner of p_map
p_map <- p_map +
  patchwork::inset_element(p_globe,
                           left = 0.58, bottom = 0.52, right = 1.0, top = 1.0,
                           align_to = "panel", clip = FALSE)

# --- Composite figures -------------------------------------------------------
sp_title <- sprintf("Historical exposure of %s", gsub("_", " ", SPECIES))
tag_theme <- theme(plot.title = element_text(face = "italic"),
                   plot.tag   = element_text(face = "bold", size = 11),
                   legend.position  = "bottom",
                   legend.direction = "horizontal")

# Standalone: variable exposure bar chart
p_combined_vars <- p_vars +
  patchwork::plot_annotation(
    title      = sp_title,
    tag_levels = list(c("a")),
    theme      = tag_theme
  )

# Standalone: heatmap + density (y-axes aligned via shared_temp_lim)
p_combined_ts <- patchwork::wrap_plots(p, p_dens,
    ncol   = 2,
    widths = c(5, 2)
  ) +
  patchwork::plot_layout(axes = "collect_x", guides = "collect") +
  patchwork::plot_annotation(
    title      = sp_title,
    tag_levels = list(c("b", "c")),
    theme      = tag_theme
  )

# Standalone: map with globe inset
p_combined_map <- p_map +
  patchwork::plot_annotation(
    title      = sp_title,
    tag_levels = list(c("d")),
    theme      = tag_theme
  )

# Full 2×2 composite: (a) p_vars, (b) p_map, (c) p heatmap, (d) p_dens
p_combined <- patchwork::wrap_plots(
    p_vars, p_map, p, p_dens,
    ncol    = 2,
    widths  = c(5, 2),
    heights = c(1.2, 3)
  ) +
  patchwork::plot_layout(axes = "collect_x") +
  patchwork::plot_annotation(
    title      = sp_title,
    tag_levels = list(c("a", "b", "", "c", "d")),
    theme      = tag_theme
  )

print(p_combined)
 
# =============================================================================
# Polar summary chart (p_polar)
# Bars = unique range cells exposed by each variable in 2025 only,
# expressed as % of total range size.
# Radius is linear 0–100%; green inner circle = "safe" zone below 25%.
# =============================================================================

all_vars_polar <- c(
  "temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo",
  "precip__12_up", "precip__12_lo", "precip__3__max_up", "precip__3__min_lo"
)

polar_raw <- tibble::tibble(var = all_vars_polar) |>
  dplyr::left_join(
    allcell_recent |>
      dplyr::filter(year == max(allcell_recent$year)) |>   # 2025 only
      dplyr::group_by(var) |>
      dplyr::summarise(n_unique = dplyr::n_distinct(cell), .groups = "drop"),
    by = "var"
  ) |>
  dplyr::mutate(
    n_unique    = tidyr::replace_na(n_unique, 0L),
    pct_range   = n_unique / range_size * 100,
    group       = dplyr::if_else(startsWith(var, "temp"), "Temperature", "Precipitation"),
    short_label = dplyr::recode(var,
      temp__12_up       = "Annual\nupper",
      temp__12_lo       = "Annual\nlower",
      temp__3__max_up   = "Seasonal\nupper",
      temp__3__min_lo   = "Seasonal\nlower",
      precip__12_up     = "Annual\nupper",
      precip__12_lo     = "Annual\nlower",
      precip__3__max_up = "Seasonal\nupper",
      precip__3__min_lo = "Seasonal\nlower"
    )
  )

# Slot order: [temp x4, GAP, precip x4, GAP] = 10 slots
# This centres temp at 12 o'clock and precip at 6 o'clock.
temp_order   <- c("temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo")
precip_order <- c("precip__3__min_lo", "precip__3__max_up", "precip__12_lo", "precip__12_up")

make_gap <- function(id) tibble::tibble(
  var = id, group = NA_character_, n_unique = NA_integer_,
  pct_range = 0, short_label = "")

polar_df <- dplyr::bind_rows(
  dplyr::filter(polar_raw, var %in% temp_order)   |> dplyr::arrange(match(var, temp_order)),
  make_gap("gap1"),
  dplyr::filter(polar_raw, var %in% precip_order) |> dplyr::arrange(match(var, precip_order)),
  make_gap("gap2")
) |>
  dplyr::mutate(
    pos      = seq_len(dplyr::n()),
    fill_col = dplyr::case_when(
      group == "Temperature"   ~ "#d73027",
      group == "Precipitation" ~ "#4575b4",
      TRUE                     ~ NA_character_
    )
  )

# In coord_polar with limits c(0.5, n_s+0.5), angle for position x:
#   theta = (x - 0.5) / n_s * 2π + start
# For temp centre (x=2.5) at 12 o'clock (theta=0): start = -(2.5-0.5)/n_s * 2π
n_s     <- nrow(polar_df)             # 10
p_start <- -(2.5 - 0.5) / n_s * 2 * pi

# Linear radius: 0–100% maps directly onto 0–max_display
ref_pcts    <- c(25, 50, 75, 100)
max_display <- 100                     # radius units = percentage points
inner_r     <- 25                      # green safe-space circle at 25%
label_r     <- max_display * 1.05     # radius for short variable labels

polar_df <- polar_df |> dplyr::mutate(r_pos = pct_range)

# Group header labels exactly at temp centre (x=2.5) and precip centre (x=7.5)
grp_labels <- tibble::tibble(
  x     = c(2.5, 7.5),
  y     = label_r * 1.08,
  label = c("TEMPERATURE", "PRECIPITATION"),
  col   = c("#d73027", "#4575b4")
)

p_polar <- ggplot(polar_df, aes(x = pos, y = r_pos)) +
  # Temperature segment background (red, 90% transparent)
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = max_display,
           fill = "#d73027", colour = NA, alpha = 0.10) +
  # Precipitation segment background (blue, 90% transparent)
  annotate("rect", xmin = 5.5, xmax = 9.5, ymin = 0, ymax = max_display,
           fill = "#4575b4", colour = NA, alpha = 0.10) +
  # Green safe-space circle in the centre (≤ 25%)
  annotate("rect", xmin = 0.5, xmax = n_s + 0.5, ymin = 0, ymax = inner_r,
           fill = "#2a9d4e", colour = NA, alpha = 0.2) +
  # Linearly-spaced reference circles
  geom_hline(yintercept = ref_pcts,
             colour = "grey80", linewidth = 0.35, linetype = "dotted") +
  # Reference circle labels at gap slot position (gap1, x = 5)
  geom_text(
    data = tibble::tibble(x = rep(5.0, length(ref_pcts)),
                          y = ref_pcts,
                          label = paste0(ref_pcts, "%")),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, size = 5, colour = "grey55", hjust = 0.5
  ) +
  # Connector lines at bar boundaries (between/around bars, not centred on them)
  geom_segment(
    data = tibble::tibble(
      xb = c(0.5, 1.5, 2.5, 3.5, 4.5,   # temp group boundaries
             5.5, 6.5, 7.5, 8.5, 9.5)    # precip group boundaries
    ),
    aes(x = xb, xend = xb, y = 0, yend = label_r * 0.88),
    colour = "grey50", linewidth = 0.35, inherit.aes = FALSE
  ) +
  # Bars
  geom_col(aes(fill = fill_col), width = 0.82,
           colour = "white", linewidth = 0.25, na.rm = TRUE) +
  # Value labels just above bars (non-zero only)
  geom_text(
    data = dplyr::filter(polar_df, pct_range > 0),
    aes(label = sprintf("%.1f%%", pct_range), y = r_pos + max_display * 0.04),
    size = 5, colour = "grey20"
  ) +
  # Short variable labels at fixed outer radius
  geom_text(
    data = dplyr::filter(polar_df, !is.na(group)),
    aes(label = short_label, y = label_r),
    size = 5, lineheight = 0.85, colour = "grey30"
  ) +
  # Group header labels
  geom_text(
    data = grp_labels,
    aes(x = x, y = y, label = label, colour = col),
    inherit.aes = FALSE, fontface = "bold", size = 5
  ) +
  coord_polar(start = p_start) +
  scale_fill_identity(na.value = NA) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0.5, n_s + 0.5), breaks = NULL) +
  scale_y_continuous(limits = c(0, label_r * 1.22), expand = c(0, 0)) +
  labs(
    title    = sprintf("%s — %% of range exposed by variable in %d",
                       gsub("_", " ", SPECIES),
                       max(allcell_recent$year)),
    subtitle = "Linear radius 0–100%; green zone = < 25%"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 15, colour = "grey30"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey50"),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1)
  )

print(p_polar)

###


