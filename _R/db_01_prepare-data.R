# =============================================================================
# 01_prepare-data.R
#
# Build script – run once during GitHub Actions (or locally) to:
#   1. Download the raw data release assets
#   2. Compute per-species exposure summaries and trend tables
#   3. Save a pre-computed cache that the Quarto dashboard reads at render time
#
# This script is intentionally transparent and step-by-step.  Each section
# is clearly labelled and outputs a log message so progress is visible in the
# GitHub Actions console.
#
# Usage (locally):
#   Rscript _R/01_prepare-data.R
#
# Environment variables (all optional; defaults listed):
#   DATA_DIR   – path where downloaded .qs data file lives (default: data/)
#   CACHE_DIR  – path to write pre-computed cache files   (default: data/)
# =============================================================================


# =============================================================================
# 0.  Load packages
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)       # fast binary serialisation (.qs files)
  library(dplyr)
  library(tidyr)
  library(purrr)     # map() for iterating over species in a functional style
})

# Source the helper functions defined in this project
source("_R/timeseries-utils.R")


# =============================================================================
# 1.  Resolve paths
# =============================================================================

data_dir  <- Sys.getenv("DATA_DIR",  "data")
cache_dir <- Sys.getenv("CACHE_DIR", "data")

# Ensure output directory exists
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Expected input file: downloaded from the GitHub release by the workflow
exposure_path <- file.path(data_dir, "allExpForShiny.qs")

stopifnot(
  "allExpForShiny.qs not found – run the download step first, or set DATA_DIR" =
    file.exists(exposure_path)
)

cat(sprintf("[prepare-data] Reading exposure data from: %s\n", exposure_path))


# =============================================================================
# 2.  Load and validate the exposure data
# =============================================================================

exposure_raw <- qs2::qs_read(exposure_path)

cat(sprintf(
  "[prepare-data] Loaded %s rows × %s columns\n",
  format(nrow(exposure_raw), big.mark = ","),
  ncol(exposure_raw)
))

# Keep only the columns the dashboard needs so the cache stays small
exposure_df <- exposure_raw |>
  dplyr::select(
    spName, cell, var, year,
    dplyr::any_of(c("propExposed", "orderName", "familyName", "redlistCategory", "group"))
  ) |>
  # Enforce correct types up-front so downstream functions don't have to guess
  dplyr::mutate(
    spName      = as.character(spName),
    var         = as.character(var),
    year        = as.integer(year),
    propExposed = as.numeric(propExposed)
  )

# Log the year range so we can confirm historical-only coverage
year_range <- range(exposure_df$year, na.rm = TRUE)
cat(sprintf(
  "[prepare-data] Year range in data: %d – %d\n",
  year_range[1], year_range[2]
))
if (year_range[2] > 2025) {
  warning(
    "[prepare-data] Data contains years beyond 2025.  ",
    "Filtering to historical-only (≤ 2025)."
  )
  exposure_df <- exposure_df |> dplyr::filter(year <= 2025)
}

# Build a catalogue of species present in the data
all_species <- sort(unique(exposure_df$spName))

# ---------------------------------------------------------------------------
# DEV MODE: set DEV_N_SPECIES to a small integer (e.g. "10") to process only
# the first N species.  Useful for fast iteration during development without
# waiting for all species to be computed.
#
# Usage:
#   DEV_N_SPECIES=10 Rscript _R/01_prepare-data.R   # shell
#   Sys.setenv(DEV_N_SPECIES = "10"); source("_R/01_prepare-data.R")  # R
#
# Leave unset (or set to "") for a full production build.
# ---------------------------------------------------------------------------
dev_n <- suppressWarnings(as.integer(Sys.getenv("DEV_N_SPECIES", "")))
if (!is.na(dev_n) && dev_n > 0) {
  all_species <- head(all_species, dev_n)
  cat(sprintf(
    "[prepare-data] DEV MODE: limiting to first %d species (set DEV_N_SPECIES= to disable)\n",
    dev_n
  ))
}

cat(sprintf("[prepare-data] Species to process: %d\n", length(all_species)))


# =============================================================================
# 3a. Load and pre-split the per-cell exposure data
#
# AllCellExposureSpXVar.rds records which cells (and variables) were exposed
# in each year.  It is split by species once here so the per-species loop
# can do O(1) lookups instead of repeated full-table scans.
# =============================================================================

allcell_path <- file.path(data_dir, "AllCellExposureSpXVar.rds")
if (file.exists(allcell_path)) {
  cat(sprintf("[prepare-data] Reading allcell data from: %s\n", allcell_path))
  allcell_df <- readRDS(allcell_path)
  allcell_year_range <- range(allcell_df$year, na.rm = TRUE)
  allcell_vars       <- sort(unique(allcell_df$var))
  cat(sprintf(
    "[prepare-data] allcell: %s rows, years %d–%d, %d variables\n",
    format(nrow(allcell_df), big.mark = ","),
    allcell_year_range[1], allcell_year_range[2],
    length(allcell_vars)
  ))
  cat("[prepare-data] Splitting allcell data by species...\n")
  allcell_list <- split(allcell_df, allcell_df$spName)
  rm(allcell_df); gc()
  cat(sprintf("[prepare-data] Split into %d species groups.\n", length(allcell_list)))
} else {
  warning("[prepare-data] AllCellExposureSpXVar.rds not found – accumulation and heatmap data will be empty.")
  allcell_list       <- list()
  allcell_year_range <- range(exposure_df$year, na.rm = TRUE)
  allcell_vars       <- sort(unique(exposure_df$var))
}


# =============================================================================
# 3b. Compute per-species exposure summaries
#
# For each species we compute five tidy tables:
#   summary      – one row per year: n_cells_exposed, mean_prop_exposed, period
#   trend        – one row per (var × year): mean_prop_exposed per variable
#   accumulation – one row per year: annual_pct, cumulative_pct of range exposed
#   var_heatmap  – one row per (var × year): pct_exposed for heatmap display
#   polar_data   – one row per variable: cumulative unique cells exposed since
#                  POLAR_YEAR_CUTOFF, expressed as % of range
#
# purrr::map() replaces an explicit for-loop; each species returns a named
# list element that becomes an entry in the cache.
# =============================================================================

POLAR_YEAR_CUTOFF <- 2023L   # first year included in polar exposure metric

cat("[prepare-data] Computing per-species timeseries (this may take a few minutes)...\n")

timeseries_cache <- purrr::map(
  .x = all_species,
  .f = function(sp) {

    # Compute overall per-year exposure summary (from propExposed data)
    summary_tbl <- compute_exposure_summary(
      species_name  = sp,
      exposure_data = exposure_df
    )

    # Compute per-variable, per-year exposure trend (from propExposed data)
    trend_tbl <- compute_species_trend(
      species_name  = sp,
      exposure_data = exposure_df
    )

    # Pre-filtered allcell rows for this species (NULL if species not in data)
    sp_allcell <- allcell_list[[sp]]

    # Compute accumulation curve (annual new + cumulative % range exposed)
    accum_tbl <- compute_exposure_accumulation(
      sp_allcell = sp_allcell,
      year_range = allcell_year_range
    )

    # Compute variable × year heatmap data
    heatmap_tbl <- compute_var_heatmap(
      sp_allcell = sp_allcell,
      year_range = allcell_year_range,
      all_vars   = allcell_vars
    )

    # Compute polar figure data: cumulative unique cells exposed per variable
    # since POLAR_YEAR_CUTOFF, as % of total range size
    if (!is.null(sp_allcell) && nrow(sp_allcell) > 0) {
      range_size <- sp_allcell$rangeSize[1]
      sp_recent  <- dplyr::filter(sp_allcell, year >= POLAR_YEAR_CUTOFF)
      polar_tbl  <- tibble::tibble(
        var = c("temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo",
                "precip__12_up", "precip__12_lo", "precip__3__max_up", "precip__3__min_lo")
      ) |>
        dplyr::left_join(
          sp_recent |>
            dplyr::group_by(var) |>
            dplyr::summarise(n_unique = dplyr::n_distinct(cell), .groups = "drop"),
          by = "var"
        ) |>
        dplyr::mutate(
          n_unique  = tidyr::replace_na(n_unique, 0L),
          pct_range = n_unique / range_size * 100
        )
    } else {
      polar_tbl <- tibble::tibble(
        var       = c("temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo",
                      "precip__12_up", "precip__12_lo", "precip__3__max_up", "precip__3__min_lo"),
        n_unique  = 0L,
        pct_range = 0
      )
    }

    list(
      summary      = summary_tbl,
      trend        = trend_tbl,
      accumulation = accum_tbl,
      var_heatmap  = heatmap_tbl,
      polar_data   = polar_tbl
    )
  }
) |>
  # Name the list elements by species name for O(1) lookup in the dashboard
  purrr::set_names(all_species)

cat(sprintf(
  "[prepare-data] Timeseries computed for %d species.\n",
  length(timeseries_cache)
))


# =============================================================================
# 4.  Build a species metadata table
#
# A lightweight lookup used by the species selector (name, order, family,
# IUCN status) without loading the full exposure_df into the dashboard.
# =============================================================================

cat("[prepare-data] Building species metadata table...\n")

species_metadata <- exposure_df |>
  dplyr::select(
    spName,
    dplyr::any_of(c("orderName", "familyName", "redlistCategory"))
  ) |>
  # One row per species (take the first non-NA value for each attribute)
  dplyr::group_by(spName) |>
  dplyr::summarise(
    dplyr::across(dplyr::everything(), ~ dplyr::first(na.omit(.x))),
    .groups = "drop"
  ) |>
  dplyr::arrange(spName)

cat(sprintf(
  "[prepare-data] Species metadata: %d rows × %d columns\n",
  nrow(species_metadata), ncol(species_metadata)
))


# =============================================================================
# 5.  Build the year-indexed exposure table
#
# The hotspot and species-level exposure tables in the dashboard filter by
# year.  Pre-splitting rows by year (as the Shiny app did) keeps the dashboard
# fast even for millions of rows.
# =============================================================================

cat("[prepare-data] Building year index...\n")

year_index <- split(seq_len(nrow(exposure_df)), exposure_df$year)

cat(sprintf(
  "[prepare-data] Year index built: %d years (%d – %d)\n",
  length(year_index),
  min(as.integer(names(year_index))),
  max(as.integer(names(year_index)))
))


# =============================================================================
# 6.  Build pre-aggregated trend tables
#
# These tables power the hotspot cell-trend chart and the per-species trend
# chart.  Pre-aggregating here avoids summarising millions of rows in the
# dashboard render.
# =============================================================================

cat("[prepare-data] Building cell-level trend table (species count per cell × year)...\n")

cell_trend_df <- exposure_df |>
  dplyr::distinct(spName, cell, year) |>   # one row per unique (sp, cell, year)
  dplyr::count(cell, year, name = "n_species_exposed")

cat("[prepare-data] Building species-level trend table (mean propExposed per sp × var × year)...\n")

sp_trend_df <- exposure_df |>
  dplyr::group_by(spName, var, year) |>
  dplyr::summarise(
    mean_prop_exposed = mean(propExposed, na.rm = TRUE),
    .groups = "drop"
  )


# =============================================================================
# 7.  Save all outputs to data/
# =============================================================================

cat("[prepare-data] Saving cache files...\n")

# Main timeseries cache – keyed by species name
cache_path <- file.path(cache_dir, "species-timeseries-cache.rds")
saveRDS(timeseries_cache, cache_path)
cat(sprintf("[prepare-data] Saved: %s  (%.1f MB)\n",
            cache_path, file.size(cache_path) / 1e6))

# Species metadata lookup
meta_path <- file.path(cache_dir, "species-metadata.rds")
saveRDS(species_metadata, meta_path)
cat(sprintf("[prepare-data] Saved: %s\n", meta_path))

# Slimmed exposure data + year index (for hotspot and species map tabs)
exposure_path_out <- file.path(cache_dir, "exposure-dashboard.rds")
saveRDS(
  list(
    exposure_df  = exposure_df,
    year_index   = year_index,
    cell_trend   = cell_trend_df,
    sp_trend     = sp_trend_df,
    years_avail  = sort(unique(exposure_df$year)),
    vars_avail   = sort(unique(exposure_df$var)),
    species_avail = all_species
  ),
  exposure_path_out
)
cat(sprintf("[prepare-data] Saved: %s  (%.1f MB)\n",
            exposure_path_out, file.size(exposure_path_out) / 1e6))

# Manifest: records what was computed and when (useful for debugging CI runs)
manifest <- list(
  computed_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
  n_species       = length(all_species),
  year_min        = year_range[1],
  year_max        = min(year_range[2], 2025L),
  source_file     = exposure_path,
  cache_files     = c(cache_path, meta_path, exposure_path_out)
)
manifest_path <- file.path(cache_dir, "build-manifest.rds")
saveRDS(manifest, manifest_path)
cat(sprintf("[prepare-data] Saved manifest: %s\n", manifest_path))


# =============================================================================
# 8.  Render polar figures as base64 PNG images
#
# Builds a planetary-boundaries-style polar chart for every species using
# the polar_data stored in timeseries_cache, then encodes each figure as a
# base64 data URI and saves to species-polar-cache.rds.
#
# This pre-rendering step keeps the dashboard render time fast.
# Requires: ggplot2, base64enc (both standard dashboard dependencies).
# =============================================================================

cat("[prepare-data] Rendering polar figures for all species...\n")

suppressPackageStartupMessages({
  library(ggplot2)
})
source("_R/timeseries-plotting.R")

polar_year_range <- c(POLAR_YEAR_CUTOFF, max(allcell_year_range))

species_figs_dir <- file.path(cache_dir, "species_figs")
dir.create(species_figs_dir, showWarnings = FALSE, recursive = TRUE)

n_polar_written <- 0L
for (sp_name in names(timeseries_cache)) {
  out_path <- file.path(species_figs_dir, paste0(sp_name, "_polar.png"))
  if (!file.exists(out_path)) {
    tryCatch({
      p <- build_polar_figure(
        polar_data = timeseries_cache[[sp_name]]$polar_data,
        sp_name    = sp_name,
        year_range = polar_year_range
      )
      ggplot2::ggsave(out_path, plot = p, width = 8, height = 8, dpi = 240)
      n_polar_written <- n_polar_written + 1L
    }, error = function(e) {
      warning(sprintf("[prepare-data] Polar figure failed for %s: %s", sp_name, e$message))
    })
  }
}
cat(sprintf("[prepare-data] Polar PNGs written to: %s (%d new)\n",
            species_figs_dir, n_polar_written))


# =============================================================================
# 9.  Build per-species range map JSON files
#
# For each species: cells exposed in the most recent year, with count of
# distinct variables exposed per cell.  Format: [[cellIdx, nVars], ...].
# OJS fetches this and renders a reactive Leaflet map; coordinates are computed
# from cellIdx using the raster geometry already passed to OJS.
# =============================================================================

cat("[prepare-data] Building species range map JSON files...\n")

map_dir <- file.path(cache_dir, "species_maps")
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

n_map_written <- 0L
for (sp in all_species) {
  out_path <- file.path(map_dir, paste0(sp, "_map.json"))
  if (!file.exists(out_path)) {
    sp_allcell <- allcell_list[[sp]]
    if (!is.null(sp_allcell) && nrow(sp_allcell) > 0) {
      recent_yr  <- max(sp_allcell$year)
      cell_vars  <- sp_allcell |>
        dplyr::filter(year == recent_yr) |>
        dplyr::group_by(cell) |>
        dplyr::summarise(n_vars = dplyr::n_distinct(var), .groups = "drop")
      # Array of [cellIdx, nVars] pairs
      jsonlite::write_json(
        lapply(seq_len(nrow(cell_vars)), function(i)
          list(as.integer(cell_vars$cell[i]), as.integer(cell_vars$n_vars[i]))),
        out_path, auto_unbox = TRUE
      )
      n_map_written <- n_map_written + 1L
    }
  }
}
cat(sprintf("[prepare-data] Range map JSONs written: %d (to %s)\n",
            n_map_written, map_dir))


# =============================================================================
# 10.  Done
# =============================================================================

cat(sprintf(
  "\n[prepare-data] Complete.\n  Species: %d\n  Years:   %d \u2013 %d\n  Output:  %s\n",
  length(all_species),
  manifest$year_min,
  manifest$year_max,
  cache_dir
))
