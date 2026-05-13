# =============================================================================
# timeseries-plotting.R
#
# Functions to produce a single combined ggplot2 figure showing historical
# species climate exposure, designed to be converted to an interactive
# plotly widget via plotly::ggplotly().
#
# The figure has two stacked panels for each species:
#   Top panel    – Mean proportion of range exposed per year (line + ribbon)
#   Bottom panel – Number of unique grid cells exposed per year (bar chart)
#
# Both panels share an x-axis (year) and are colour-coded by time period
# (Historical baseline vs. Recent observed).
#
# The output of plot_species_timeseries_combined() is a standard ggplot2 object.
# The caller wraps it in plotly::ggplotly() for interactivity.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)


# Colour palette for time periods – distinct, colourblind-friendly
PERIOD_COLOURS <- c(
  "Historical baseline" = "#4575b4",   # blue
  "Recent observed"     = "#d73027"    # red
)

# Colour palette per variable – reds for temperature, blues for precipitation
VAR_COLOURS <- c(
  temp__12_up       = "#d73027",
  temp__12_lo       = "#fc8d59",
  temp__3__max_up   = "#fdae61",
  temp__3__min_lo   = "#f46d43",
  precip__12_up     = "#313695",
  precip__12_lo     = "#4575b4",
  precip__3__max_up = "#74add1",
  precip__3__min_lo = "#abd9e9"
)

# Human-readable variable labels (matches VAR_LABELS in app.R)
VAR_LABELS <- c(
  temp__12_up       = "Temp: annual upper extreme",
  temp__12_lo       = "Temp: annual lower extreme",
  temp__3__max_up   = "Temp: seasonal max upper",
  temp__3__min_lo   = "Temp: seasonal min lower",
  precip__12_up     = "Precip: annual upper extreme",
  precip__12_lo     = "Precip: annual lower extreme",
  precip__3__max_up = "Precip: seasonal max upper",
  precip__3__min_lo = "Precip: seasonal min lower"
)

# Helper: translate raw variable codes to display labels
label_var <- function(var_code) {
  lab <- VAR_LABELS[as.character(var_code)]
  dplyr::if_else(is.na(lab), as.character(var_code), lab)
}


# -----------------------------------------------------------------------------
# plot_species_timeseries_combined()
#
# Produces a single combined ggplot2 figure with two vertically stacked panels
# for a focal species, using only historical data.
#
# Arguments:
#   exposure_summary – data frame from compute_exposure_summary(); columns:
#                      year, n_cells_exposed, mean_prop_exposed, period
#   species_trend    – data frame from compute_species_trend(); columns:
#                      var, year, mean_prop_exposed, period
#   species_name     – character; used in the figure title
#
# Returns:
#   A ggplot2 object suitable for plotly::ggplotly().
#   Facet rows: "% range exposed" (top) and "Cells exposed" (bottom).
# -----------------------------------------------------------------------------
plot_species_timeseries_combined <- function(exposure_summary,
                                             species_trend,
                                             species_name) {

  # --- 1. Prepare overall exposure data for the top panel -------------------
  # Convert mean_prop_exposed to a percentage so the y-axis is intuitive.
  overall <- exposure_summary |>
    dplyr::mutate(
      pct_exposed = mean_prop_exposed * 100,
      # Facet label that will appear as the row heading in the combined figure
      panel = "% range exposed (mean across variables)"
    )

  # --- 2. Prepare per-variable trend data for the bottom panel --------------
  # Translate raw variable codes to readable labels for display in the legend.
  by_var <- species_trend |>
    dplyr::mutate(
      var_label = label_var(var),
      pct_exposed = mean_prop_exposed * 100,
      panel = "% range exposed (by variable)"
    )

  # --- 3. Build the combined figure using facet_grid ------------------------
  # We use two separate geom layers with the same x/y aesthetics so that
  # ggplotly can render both panels in one interactive widget.
  #
  # Facet rows are driven by the "panel" column present in both data frames.

  ggplot() +

    # Top panel: overall mean proportion of range exposed through time.
    # geom_area fills the region under the curve so total exposure through
    # time is immediately visible; geom_line adds a clear trend line.
    geom_area(
      data    = overall,
      mapping = aes(x = year, y = pct_exposed, fill = period,
                    text = paste0(
                      "Year: ", year,
                      "\nPeriod: ", period,
                      "\n% range exposed: ", round(pct_exposed, 2), "%"
                    )),
      alpha = 0.4
    ) +
    geom_line(
      data    = overall,
      mapping = aes(x = year, y = pct_exposed, colour = period),
      linewidth = 0.8
    ) +

    # Bottom panel: per-variable mean proportion of range exposed.
    # Separate lines per variable allow users to see which climate driver
    # is responsible for exposure in any given year.
    geom_line(
      data    = by_var,
      mapping = aes(x = year, y = pct_exposed, colour = var_label,
                    group = var_label,
                    text = paste0(
                      "Variable: ", var_label,
                      "\nYear: ", year,
                      "\n% range exposed: ", round(pct_exposed, 2), "%"
                    )),
      linewidth = 0.6,
      alpha = 0.85
    ) +

    # --- Split into two rows by the "panel" column -------------------------
    facet_grid(
      rows = vars(panel),
      scales = "free_y"   # each panel can have its own y-axis range
    ) +

    # --- Axes and labels ---------------------------------------------------
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8)
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%")
    ) +

    # Period colours for the overall panel; variable colours handled automatically
    scale_colour_manual(
      values = PERIOD_COLOURS,
      na.value = "grey60"
    ) +
    scale_fill_manual(
      values = PERIOD_COLOURS,
      na.value = "grey80"
    ) +

    # --- Labels ------------------------------------------------------------
    labs(
      title   = paste0("Historical climate exposure: ", species_name),
      x       = "Year",
      y       = "% of range exposed",
      colour  = NULL,
      fill    = NULL,
      caption = paste0(
        "Data: allExpForShiny  |  ",
        "Historical baseline: up to 2010  |  Recent observed: 2011\u20132025"
      )
    ) +

    # --- Theme: clean scientific style ------------------------------------
    theme_bw(base_size = 12) +
    theme(
      # Move legend below plot so it does not compress the panels
      legend.position  = "bottom",
      legend.direction = "horizontal",
      # Facet strip labels are the panel headings
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text       = element_text(size = 10, face = "bold"),
      # Minimal grid lines – keep y-grid for readability, remove x-grid
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      # Ensure titles are readable at dashboard size
      plot.title   = element_text(size = 13, face = "bold"),
      plot.caption = element_text(size = 8, colour = "grey50")
    )
}


# -----------------------------------------------------------------------------
# plot_cumulative_exposure()
#
# For the "Recent observed" period (year > 2010), plots:
#   Grey area  – cumulative % of range exposed (running sum of annual mean % exposed)
#   Coloured lines – per-variable % of range exposed each year
#
# Arguments:
#   exposure_summary – from compute_exposure_summary(); used for cumulative area
#   species_trend    – from compute_species_trend(); used for per-variable lines
#   species_name     – character; figure title
#
# Per-variable trace names are encoded as "label||var_code" so OJS can parse
# the var_code for checkbox filtering without losing the display label.
# -----------------------------------------------------------------------------
plot_cumulative_exposure <- function(exposure_summary, species_trend = NULL, species_name) {

  recent_summary <- exposure_summary |>
    dplyr::filter(period == "Recent observed") |>
    dplyr::arrange(year) |>
    dplyr::mutate(cumulative_pct = pmin(cumsum(mean_prop_exposed * 100), 100))

  p <- plotly::plot_ly() |>
    # Grey shaded area: cumulative % of range exposed (left axis)
    plotly::add_trace(
      data          = recent_summary,
      x             = ~year,
      y             = ~cumulative_pct,
      type          = "scatter",
      mode          = "lines",
      fill          = "tozeroy",
      fillcolor     = "rgba(150,150,150,0.35)",
      line          = list(color = "rgba(100,100,100,0.6)", width = 1),
      name          = "Cumulative % exposed",
      yaxis         = "y",
      hovertemplate = "Year: %{x}<br>Cumulative % exposed: %{y:.1f}%<extra></extra>"
    )

  # Per-variable lines (right axis, % range exposed) – filtered by OJS checkboxes
  if (!is.null(species_trend)) {
    recent_trend <- species_trend |>
      dplyr::filter(year > 2010) |>
      dplyr::arrange(year)

    vars <- unique(recent_trend$var)

    for (v in vars) {
      vdata  <- dplyr::filter(recent_trend, var == v)
      label  <- VAR_LABELS[v]
      if (is.na(label)) label <- v
      colour <- VAR_COLOURS[v]
      if (is.na(colour)) colour <- "#888888"

      # Encode var_code in trace name so OJS can parse it for checkbox filtering
      trace_name <- paste0(label, "||", v)

      p <- plotly::add_trace(p,
        data          = vdata,
        x             = ~year,
        y             = ~I(mean_prop_exposed * 100),
        type          = "scatter",
        mode          = "lines+markers",
        line          = list(color = colour, width = 1.8),
        marker        = list(color = colour, size = 5),
        name          = trace_name,
        yaxis         = "y2",
        hovertemplate = paste0(label, "<br>Year: %{x}<br>% exposed: %{y:.1f}%<extra></extra>")
      )
    }
  }

  plotly::layout(p,
    title  = list(
      text = paste0("Recent exposure accumulation: ", species_name),
      font = list(size = 14, family = "sans-serif")
    ),
    xaxis  = list(title = "Year"),
    yaxis  = list(title = "Cumulative % of range exposed", side = "left", ticksuffix = "%"),
    yaxis2 = list(
      title      = "% range exposed",
      side       = "right",
      overlaying = "y",
      showgrid   = FALSE,
      ticksuffix = "%"
    ),
    legend    = list(orientation = "h", y = -0.25),
    hovermode = "x unified"
  )
}


# -----------------------------------------------------------------------------
# plot_species_exposure_figure()
#
# Produces a two-panel interactive plotly figure for a single species:
#
#   Top panel    – Heatmap of % range exposed per climate variable × year.
#                  Each row is one climate variable; colour encodes the
#                  percentage of range cells exposed that year.
#                  Colour scale: light yellow (0%) → orange → red →
#                  dark purple (100%), matching the original figure style.
#
#   Bottom panel – Cumulative exposure accumulation.
#                  Grey filled area: % of range ever exposed (non-decreasing).
#                  Red bars: % of range cells first exposed each year.
#
# Arguments:
#   var_heatmap   – from compute_var_heatmap(); columns: var, year, pct_exposed
#   accumulation  – from compute_exposure_accumulation(); columns: year,
#                   annual_pct, cumulative_pct
#   species_name  – character; used in figure title
#
# Returns a native plotly object (subplot of two panels).
# -----------------------------------------------------------------------------
plot_species_exposure_figure <- function(var_heatmap, accumulation, species_name) {

  # Variable display order: temperature upper/lower, then precipitation
  var_order <- c(
    "temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo",
    "precip__12_up", "precip__12_lo", "precip__3__max_up", "precip__3__min_lo"
  )
  var_order <- var_order[var_order %in% unique(var_heatmap$var)]

  # Pivot to wide matrix: rows = variables, columns = years
  mat_wide <- var_heatmap |>
    dplyr::mutate(var = factor(var, levels = var_order)) |>
    dplyr::filter(!is.na(var)) |>
    dplyr::arrange(var, year) |>
    tidyr::pivot_wider(
      id_cols      = var,
      names_from   = year,
      values_from  = pct_exposed,
      values_fill  = 0
    )

  # Shorten labels for the y-axis
  y_labs <- VAR_LABELS[as.character(mat_wide$var)]
  y_labs[is.na(y_labs)] <- as.character(mat_wide$var)[is.na(y_labs)]
  y_labs <- gsub("Temperature: ", "Temp: ", gsub("Precipitation: ", "Precip: ", y_labs))

  mat_vals <- as.matrix(mat_wide[, -1, drop = FALSE])
  rownames(mat_vals) <- y_labs
  years_x  <- as.integer(colnames(mat_wide)[-1])

  # Custom colour scale: light yellow → orange → red → dark purple
  cscale <- list(
    list(0,   "#ffffcc"),
    list(0.2, "#fd8d3c"),
    list(0.4, "#e31a1c"),
    list(0.7, "#800026"),
    list(1.0, "#4d004b")
  )

  # ── Top panel: variable × year heatmap ─────────────────────────────────────
  p_top <- plotly::plot_ly(
    z             = mat_vals,
    x             = years_x,
    y             = rownames(mat_vals),
    type          = "heatmap",
    colorscale    = cscale,
    zmin          = 0,
    zmax          = 100,
    colorbar      = list(
      title  = "% exposed",
      ticksuffix = "%",
      len    = 0.5,
      y      = 0.78
    ),
    hovertemplate = "Year: %{x}<br>%{y}<br>% exposed: %{z:.1f}%<extra></extra>"
  )

  # ── Bottom panel: cumulative area + annual bars ─────────────────────────────
  p_bottom <- plotly::plot_ly() |>

    # Grey filled area – cumulative % of range ever exposed
    plotly::add_trace(
      data          = accumulation,
      x             = ~year,
      y             = ~cumulative_pct,
      type          = "scatter",
      mode          = "lines",
      fill          = "tozeroy",
      fillcolor     = "rgba(150,150,150,0.5)",
      line          = list(color = "rgba(80,80,80,0.8)", width = 1),
      name          = "Cumulative % exposed",
      hovertemplate = "Year: %{x}<br>Cumulative: %{y:.1f}%<extra></extra>"
    ) |>

    # Red bars – % of range cells newly exposed each year
    plotly::add_trace(
      data          = dplyr::filter(accumulation, annual_pct > 0),
      x             = ~year,
      y             = ~annual_pct,
      type          = "bar",
      marker        = list(color = "#d73027"),
      name          = "Annual % newly exposed",
      hovertemplate = "Year: %{x}<br>Annual new: %{y:.1f}%<extra></extra>"
    )

  # ── Combine into a 2-row subplot ───────────────────────────────────────────
  plotly::subplot(
    p_top, p_bottom,
    nrows   = 2,
    shareX  = TRUE,
    titleY  = TRUE,
    heights = c(0.42, 0.58)
  ) |>
    plotly::layout(
      title     = list(
        text = paste0("Climate exposure: ", gsub("_", " ", species_name)),
        font = list(size = 13, family = "sans-serif")
      ),
      xaxis  = list(title = "Year"),
      yaxis  = list(title = "Climate variable"),
      yaxis2 = list(
        title      = "% range exposed",
        ticksuffix = "%",
        range      = c(0, 100)
      ),
      barmode   = "overlay",
      hovermode = "x unified",
      legend    = list(orientation = "h", y = -0.12)
    )
}


# =============================================================================
# fig_to_b64()
#
# Renders a ggplot2 object to a PNG and returns a base64 data URI string
# suitable for embedding in HTML <img src="..."> tags.
#
# Arguments:
#   p      – ggplot2 object
#   width  – figure width in inches  (default 5)
#   height – figure height in inches (default 5)
#   dpi    – resolution (default 120)
#
# Returns a character string: "data:image/png;base64,<encoded data>"
# =============================================================================
fig_to_b64 <- function(p, width = 5, height = 5, dpi = 120) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  ggplot2::ggsave(tmp, plot = p, width = width, height = height,
                  dpi = dpi, bg = "white")
  raw_bytes <- readBin(tmp, what = "raw", n = file.info(tmp)$size)
  b64 <- base64enc::base64encode(raw_bytes)
  paste0("data:image/png;base64,", b64)
}


# =============================================================================
# build_polar_figure()
#
# Builds a planetary-boundaries-style polar bar chart showing the cumulative
# percentage of a species' range exposed to each climate variable.
#
# Arguments:
#   polar_data  – data frame with columns: var (character), pct_range (numeric)
#                 One row per variable (8 rows total).  Produced by the build
#                 step in 01_prepare-data.R.
#   sp_name     – species name string (used in plot title)
#   year_range  – integer vector of length 2: c(min_year, max_year) for subtitle
#
# Returns a ggplot2 object.
# =============================================================================
build_polar_figure <- function(polar_data, sp_name, year_range) {

  all_vars_polar <- c(
    "temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo",
    "precip__12_up", "precip__12_lo", "precip__3__max_up", "precip__3__min_lo"
  )

  short_labels <- c(
    temp__12_up       = "Annual\nupper",  temp__12_lo       = "Annual\nlower",
    temp__3__max_up   = "Seasonal\nupper", temp__3__min_lo  = "Seasonal\nlower",
    precip__12_up     = "Annual\nupper",  precip__12_lo     = "Annual\nlower",
    precip__3__max_up = "Seasonal\nupper", precip__3__min_lo = "Seasonal\nlower"
  )

  polar_raw <- tibble::tibble(var = all_vars_polar) |>
    dplyr::left_join(
      dplyr::select(polar_data, var, pct_range),
      by = "var"
    ) |>
    dplyr::mutate(
      pct_range   = tidyr::replace_na(pct_range, 0),
      group       = dplyr::if_else(startsWith(var, "temp"), "Temperature", "Precipitation"),
      short_label = short_labels[var]
    )

  temp_order   <- c("temp__12_up", "temp__12_lo", "temp__3__max_up", "temp__3__min_lo")
  precip_order <- c("precip__3__min_lo", "precip__3__max_up", "precip__12_lo", "precip__12_up")

  make_gap <- function(id) tibble::tibble(
    var = id, group = NA_character_, pct_range = 0, short_label = "")

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

  n_s         <- nrow(polar_df)
  p_start     <- -(2.5 - 0.5) / n_s * 2 * pi
  ref_pcts    <- c(0.1, 1, 10, 50)
  max_display <- 20
  to_r        <- function(pct) log1p(pct * 10) / log1p(100 * 10) * max_display
  inner_r     <- to_r(0.2)
  label_r     <- max_display * 1.05

  polar_df <- polar_df |> dplyr::mutate(r_pos = to_r(pct_range))

  grp_labels <- tibble::tibble(
    x     = c(2.5, 7.5),
    y     = label_r * 1.08,
    label = c("TEMPERATURE", "PRECIPITATION"),
    col   = c("#d73027", "#4575b4")
  )

  ggplot2::ggplot(polar_df, ggplot2::aes(x = pos, y = r_pos)) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = max_display,
                      fill = "#d73027", colour = NA, alpha = 0.10) +
    ggplot2::annotate("rect", xmin = 5.5, xmax = 9.5, ymin = 0, ymax = max_display,
                      fill = "#4575b4", colour = NA, alpha = 0.10) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = n_s + 0.5, ymin = 0, ymax = inner_r,
                      fill = "#2a9d4e", colour = NA, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = to_r(ref_pcts),
                        colour = "grey80", linewidth = 0.35, linetype = "dotted") +
    ggplot2::geom_text(
      data = tibble::tibble(x = rep(5.0, length(ref_pcts)), y = to_r(ref_pcts),
                            label = paste0(ref_pcts, "%")),
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, size = 2.5, colour = "grey55", hjust = 0.5
    ) +
    ggplot2::geom_segment(
      data = tibble::tibble(xb = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)),
      ggplot2::aes(x = xb, xend = xb, y = 0, yend = label_r * 0.88),
      colour = "grey50", linewidth = 0.35, inherit.aes = FALSE
    ) +
    ggplot2::geom_col(ggplot2::aes(fill = fill_col), width = 0.82,
                      colour = "white", linewidth = 0.25, na.rm = TRUE) +
    ggplot2::geom_text(
      data = dplyr::filter(polar_df, pct_range > 0),
      ggplot2::aes(label = sprintf("%.1f%%", pct_range), y = r_pos + max_display * 0.04),
      size = 2.5, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(polar_df, !is.na(group)),
      ggplot2::aes(label = short_label, y = label_r),
      size = 2.5, lineheight = 0.85, colour = "grey30"
    ) +
    ggplot2::geom_text(
      data = grp_labels,
      ggplot2::aes(x = x, y = y, label = label, colour = col),
      inherit.aes = FALSE, fontface = "bold", size = 3.5
    ) +
    ggplot2::coord_polar(start = p_start) +
    ggplot2::scale_fill_identity(na.value = NA) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(limits = c(0.5, n_s + 0.5), breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, label_r * 1.22), expand = c(0, 0)) +
    ggplot2::labs(
      title    = sprintf("%s \u2014 cumulative %% of range exposed (%d\u2013%d)",
                         gsub("_", " ", sp_name), year_range[1], year_range[2]),
      subtitle = "Radius is log-transformed"
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey30"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 7.5, colour = "grey50"),
      plot.margin   = ggplot2::margin(1, 1, 1, 1)
    )
}
