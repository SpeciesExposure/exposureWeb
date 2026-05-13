 # =============================================================================
# timeseries-utils.R
#
# Functions to compute historical species climate exposure summaries from the
# allExpForShiny dataset.  All data is strictly historical (up to 2025); there
# are no forecast or projection steps.
#
# The data model:
#   allExpForShiny has one row per (species × cell × climate variable × year).
#   "propExposed" is the fraction of that cell's area that exceeded the
#   climate threshold for that variable in that year.
#
# Functions are tidyverse-first, written for clarity over speed, and designed
# to produce tidy data frames that feed directly into ggplot2 / ggplotly.
# =============================================================================

library(dplyr)
library(tidyr)


# -----------------------------------------------------------------------------
# classify_time_period()
#
# Assigns a human-readable period label to each year based on the data
# collection context used in this study:
#   "Historical baseline"  – 1980–2010  (climatological reference period)
#   "Recent observed"      – 2011–2025  (monitoring window, fully observed)
#
# Arguments:
#   year  – integer vector of years
#
# Returns a character vector of period labels (same length as year).
# -----------------------------------------------------------------------------
classify_time_period <- function(year) {
  dplyr::case_when(
    year <= 2020 ~ "Historical baseline",
    year <= 2025 ~ "Recent observed",
    TRUE         ~ NA_character_   # safety catch; no future years expected
  )
}


# -----------------------------------------------------------------------------
# compute_exposure_summary()
#
# For a single species, computes per-year exposure summaries across all
# climate variables and grid cells present in the data.
#
# "Exposure" here means: for a given year, how many unique grid cells contained
# at least one climate threshold exceedance, and what was the average
# propExposed across those cells × variable combinations?
#
# Arguments:
#   species_name   – character; must match a value in app_df$spName
#   exposure_data  – data frame with columns: spName, cell, var, year, propExposed
#
# Returns a tidy data frame with one row per year:
#   year             – integer
#   n_cells_exposed  – number of unique grid cells exposed that year
#   mean_prop_exposed – mean propExposed across all (cell × var) rows for that year
#   period           – time period label (from classify_time_period)
# -----------------------------------------------------------------------------
compute_exposure_summary <- function(species_name, exposure_data) {

  exposure_data |>
    # Keep only rows for the focal species
    dplyr::filter(spName == species_name) |>

    # Ensure propExposed is numeric (it may be stored as character in some builds)
    dplyr::mutate(propExposed = as.numeric(propExposed)) |>

    # Summarise per year: count unique exposed cells and average exposure intensity
    dplyr::group_by(year) |>
    dplyr::summarise(
      n_cells_exposed   = dplyr::n_distinct(cell),
      mean_prop_exposed = mean(propExposed, na.rm = TRUE),
      .groups = "drop"
    ) |>

    # Fill in years with zero exposure so the timeline is continuous
    # (species may not appear in allExpForShiny for years with no exposure)
    tidyr::complete(
      year = seq(min(exposure_data$year), max(exposure_data$year)),
      fill = list(n_cells_exposed = 0L, mean_prop_exposed = 0)
    ) |>

    # Attach human-readable period labels
    dplyr::mutate(period = classify_time_period(year)) |>

    dplyr::arrange(year)
}


# -----------------------------------------------------------------------------
# compute_species_trend()
#
# For a single species, computes the per-variable, per-year mean proportion of
# range exposed.  This feeds the per-variable breakdown panel in the timeseries
# figure so users can see which climate variables are driving exposure.
#
# Arguments:
#   species_name   – character; must match a value in app_df$spName
#   exposure_data  – data frame with columns: spName, var, year, propExposed
#
# Returns a tidy data frame with one row per (var × year):
#   var              – climate variable code
#   year             – integer
#   mean_prop_exposed – mean propExposed for that species × var × year
#   period           – time period label
# -----------------------------------------------------------------------------
compute_species_trend <- function(species_name, exposure_data) {

  # Get the full year range from the dataset so we can fill gaps
  year_range <- seq(min(exposure_data$year), max(exposure_data$year))
  all_vars   <- unique(exposure_data$var)

  exposure_data |>
    dplyr::filter(spName == species_name) |>
    dplyr::mutate(propExposed = as.numeric(propExposed)) |>

    # Average propExposed across cells for each variable × year combination
    dplyr::group_by(var, year) |>
    dplyr::summarise(
      mean_prop_exposed = mean(propExposed, na.rm = TRUE),
      .groups = "drop"
    ) |>

    # Fill every (var × year) combination, putting 0 where no exposure occurred
    tidyr::complete(
      var  = all_vars,
      year = year_range,
      fill = list(mean_prop_exposed = 0)
    ) |>

    dplyr::mutate(period = classify_time_period(year)) |>
    dplyr::arrange(var, year)
}


# -----------------------------------------------------------------------------
# compute_exposure_accumulation()
#
# For a single species, computes per-year counts of *newly* exposed cells
# (first time a cell exceeded any threshold) and the running cumulative
# percentage of the range ever exposed.
#
# Arguments:
#   sp_allcell  – pre-filtered rows from AllCellExposureSpXVar for one species;
#                 columns: cell, year, rangeSize, var, spName, group
#   year_range  – integer vector of length 2: c(min_year, max_year)
#
# Returns a tidy data frame with one row per year:
#   year           – integer
#   annual_pct     – % of range cells first exposed in that year
#   cumulative_pct – running % of range ever exposed (capped at 100)
#   period         – time period label
# -----------------------------------------------------------------------------
compute_exposure_accumulation <- function(sp_allcell, year_range) {

  if (is.null(sp_allcell) || nrow(sp_allcell) == 0) {
    yrs <- seq(year_range[1], year_range[2])
    return(tibble::tibble(
      year           = yrs,
      annual_pct     = 0,
      cumulative_pct = 0,
      period         = classify_time_period(yrs)
    ))
  }

  range_size <- sp_allcell$rangeSize[1]

  # First year each cell was ever exposed (any variable)
  first_exp <- sp_allcell |>
    dplyr::group_by(cell) |>
    dplyr::summarise(first_year = min(year), .groups = "drop") |>
    dplyr::count(first_year, name = "n_new_cells")

  tibble::tibble(year = seq(year_range[1], year_range[2])) |>
    dplyr::left_join(first_exp, by = c("year" = "first_year")) |>
    dplyr::mutate(
      n_new_cells    = dplyr::coalesce(n_new_cells, 0L),
      cumulative_pct = pmin(cumsum(n_new_cells) / range_size * 100, 100),
      annual_pct     = n_new_cells / range_size * 100,
      period         = classify_time_period(year)
    ) |>
    dplyr::select(year, annual_pct, cumulative_pct, period)
}


# -----------------------------------------------------------------------------
# compute_var_heatmap()
#
# For a single species, computes per-variable × per-year percentage of range
# cells exposed, for use as a heatmap in the exposure figure.
#
# Arguments:
#   sp_allcell  – pre-filtered rows from AllCellExposureSpXVar for one species
#   year_range  – integer vector of length 2: c(min_year, max_year)
#   all_vars    – character vector of all variable codes (to fill in zeros)
#
# Returns a tidy data frame with one row per (var × year):
#   var         – climate variable code
#   year        – integer
#   pct_exposed – % of range cells exposed to that variable in that year
#   period      – time period label
# -----------------------------------------------------------------------------
compute_var_heatmap <- function(sp_allcell, year_range, all_vars) {

  if (is.null(sp_allcell) || nrow(sp_allcell) == 0) {
    return(
      tidyr::expand_grid(
        var  = all_vars,
        year = seq(year_range[1], year_range[2])
      ) |>
        dplyr::mutate(pct_exposed = 0, period = classify_time_period(year))
    )
  }

  range_size <- sp_allcell$rangeSize[1]

  sp_allcell |>
    dplyr::group_by(var, year) |>
    dplyr::summarise(n_cells = dplyr::n_distinct(cell), .groups = "drop") |>
    tidyr::complete(
      var  = all_vars,
      year = seq(year_range[1], year_range[2]),
      fill = list(n_cells = 0L)
    ) |>
    dplyr::mutate(
      pct_exposed = n_cells / range_size * 100,
      period      = classify_time_period(year)
    ) |>
    dplyr::select(var, year, pct_exposed, period) |>
    dplyr::arrange(var, year)
}
