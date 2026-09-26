# Shared visual design and forest plots

#' Use the ekbMed plot theme
#'
#' Provides a compact clinical-reporting ggplot theme.
#'
#' @param base_size Base plot font size in points.
#' @param base_family Base plot font family.
#' @return A ggplot2 theme.
#' @export
#' @examples
#' theme_ekbmed(base_size = 11)
theme_ekbmed <- function(base_size = 11, base_family = "") {
  .ekb_require("ggplot2")
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(hjust = 0)
    )
}

# Internal journal-inspired colour palettes
.ekb_palette <- function(
    palette = c("lancet", "nejm", "jama"),
    n = 2
) {
  
  palette <- match.arg(palette)
  
  values <- switch(
    palette,
    
    lancet = c(
      "#00468B",
      "#ED0000",
      "#42B540",
      "#0099B4",
      "#925E9F",
      "#FDAF91",
      "#AD002A",
      "#ADB6B6",
      "#1B1919"
    ),
    
    nejm = c(
      "#BC3C29",
      "#0072B5",
      "#E18727",
      "#20854E",
      "#7876B1",
      "#6F99AD",
      "#FFDC91",
      "#EE4C97"
    ),
    
    jama = c(
      "#374E55",
      "#DF8F44",
      "#00A1D5",
      "#B24745",
      "#79AF97",
      "#6A6599",
      "#80796B"
    )
  )
  
  if (n > length(values)) {
    stop(
      sprintf(
        "Palette '%s' contains only %s colours.",
        palette,
        length(values)
      ),
      call. = FALSE
    )
  }
  
  values[seq_len(n)]
}

#' Apply the legacy ekbMed table theme
#'
#' Styles a flextable; use as_ekb_flextable for the publication-table API.
#'
#' @param ft A flextable object.
#' @param textsize Body text size in points.
#' @param headsize Header text size in points.
#' @return A flextable object.
#' @export
#' @examples
#' if (requireNamespace("flextable", quietly = TRUE)) {
#'   theme_ekbmed_flextable(flextable::flextable(data.frame(a = 1:3)))
#' }
theme_ekbmed_flextable <- function(ft, textsize = 9, headsize = 10) {
  .ekb_require("flextable")
  ft <- flextable::autofit(ft)
  ft <- flextable::fontsize(ft, size = textsize, part = "body")
  ft <- flextable::fontsize(ft, size = headsize, part = "header")
  ft <- flextable::bold(ft, bold = TRUE, part = "header")
  ft <- flextable::padding(ft, padding = 2, part = "all")
  ft <- flextable::line_spacing(ft, space = 1, part = "all")
  ft <- flextable::align(ft, align = "left", part = "header")
  ft
}

#' Scale a flextable to page width
#'
#' Scales table column widths to a specified printable width.
#'
#' @param ft A flextable object.
#' @param pgwidth Available table width in inches.
#' @return A flextable object.
#' @export
#' @examples
#' if (requireNamespace("flextable", quietly = TRUE)) {
#'   fit_flextable_to_page(flextable::flextable(data.frame(a = 1:3)))
#' }
fit_flextable_to_page <- function(ft, pgwidth = 6) {
  .ekb_require("flextable")
  ft <- flextable::autofit(ft)
  if (length(pgwidth) != 1L || !is.finite(pgwidth) || pgwidth <= 0) stop("pgwidth must be positive.")
  dims <- flextable::dim_pretty(ft)
  if (!is.null(dims$widths) && sum(dims$widths) > 0) {
    ft <- flextable::width(ft, width = dims$widths * pgwidth / sum(dims$widths))
  }
  ft
}


# ============================================================
# Forest plot helpers
# ============================================================


.ekb_forest_label <- function(x, labels = NULL) {
  
  x <- as.character(x)
  
  if (is.null(labels)) {
    return(x)
  }
  
  if (is.null(names(labels))) {
    stop(
      "labels must be a named character vector.",
      call. = FALSE
    )
  }
  
  if (x %in% names(labels)) {
    return(as.character(labels[[x]]))
  }
  
  x
}


.ekb_forest_level_label <- function(
    variable,
    level,
    level_labels = NULL
) {
  
  level <- as.character(level)
  
  if (is.null(level_labels)) {
    return(level)
  }
  
  mapping <- level_labels[[variable]]
  
  if (is.null(mapping)) {
    return(level)
  }
  
  if (is.null(names(mapping))) {
    stop(
      sprintf(
        "level_labels[['%s']] must be a named character vector.",
        variable
      ),
      call. = FALSE
    )
  }
  
  if (level %in% names(mapping)) {
    return(as.character(mapping[[level]]))
  }
  
  level
}


.ekb_forest_limits <- function(
    conf.low,
    conf.high,
    limits = NULL
) {
  
  if (!is.null(limits)) {
    
    if (
      length(limits) != 2L ||
      any(!is.finite(limits)) ||
      any(limits <= 0) ||
      limits[1] >= limits[2]
    ) {
      stop(
        "limits must contain two positive increasing values.",
        call. = FALSE
      )
    }
    
    return(as.numeric(limits))
  }
  
  
  x <- c(conf.low, conf.high)
  
  x <- x[
    is.finite(x) &
      x > 0
  ]
  
  if (!length(x)) {
    return(c(0.5, 2))
  }
  
  
  lo <- min(c(x, 1))
  hi <- max(c(x, 1))
  
  
  candidates <- c(
    0.02,
    0.05,
    0.1,
    0.2,
    0.25,
    0.5,
    1,
    2,
    4,
    5,
    10,
    20,
    50
  )
  
  
  lower_candidates <- candidates[
    candidates <= lo
  ]
  
  upper_candidates <- candidates[
    candidates >= hi
  ]
  
  
  lower <- if (length(lower_candidates)) {
    max(lower_candidates)
  } else {
    10 ^ floor(log10(lo))
  }
  
  
  upper <- if (length(upper_candidates)) {
    min(upper_candidates)
  } else {
    10 ^ ceiling(log10(hi))
  }
  
  
  # Keep HR = 1 away from the plot boundary
  if (lower >= 1) {
    lower <- 0.5
  }
  
  if (upper <= 1) {
    upper <- 2
  }
  
  
  c(lower, upper)
}


.ekb_forest_breaks <- function(
    limits,
    breaks = NULL
) {
  
  if (!is.null(breaks)) {
    
    breaks <- as.numeric(breaks)
    
    breaks <- breaks[
      is.finite(breaks) &
        breaks > 0 &
        breaks >= limits[1] &
        breaks <= limits[2]
    ]
    
    return(sort(unique(breaks)))
  }
  
  
  candidates <- c(
    0.02,
    0.05,
    0.1,
    0.2,
    0.25,
    0.5,
    1,
    2,
    4,
    5,
    10,
    20,
    50
  )
  
  out <- candidates[
    candidates >= limits[1] &
      candidates <= limits[2]
  ]
  
  
  if (!1 %in% out) {
    out <- sort(unique(c(out, 1)))
  }
  
  
  out
}


# ============================================================
# Internal publication-style forest plot
# ============================================================

.ekb_plot_forest <- function(
    data,
    show_n = FALSE,
    show_p = TRUE,
    variable_heading = "Variable",
    level_heading = "Category",
    n_heading = "N",
    hr_heading = "HR (95% CI)",
    p_heading = "P value",
    title = NULL,
    xlab = "Hazard ratio",
    limits = NULL,
    breaks = NULL,
    reference = 1,
    reference_label = "Reference",
    digits = 2,
    p_digits = 3,
    base_size = 10,
    point_size = 2.5,
    ci_linewidth = 0.55,
    reference_linewidth = 0.45,
    block_gap = 0.12,
    shade = TRUE,
    column_widths = NULL,
    show_favors = FALSE,
    favors_left = NULL,
    favors_right = NULL
) {
  
  .ekb_require("ggplot2")
  
  d <- as.data.frame(data)
  
  # ----------------------------------------------------------
  # Typography
  # ----------------------------------------------------------
  
  text_size <- base_size / 3.3
  header_size <- base_size / 3.15
  
  axis_text_size <- text_size
  axis_title_size <- header_size
  
  # ----------------------------------------------------------
  # Validate input
  # ----------------------------------------------------------
  
  required <- c(
    "block",
    "variable",
    "level",
    "estimate",
    "conf.low",
    "conf.high",
    "p.value",
    "n",
    "is_reference"
  )
  
  .ekb_assert_columns(
    d,
    required
  )
  
  if (!nrow(d)) {
    stop(
      "No data available for forest plot.",
      call. = FALSE
    )
  }
  
  
  # ----------------------------------------------------------
  # Display strings
  # ----------------------------------------------------------
  
  d$variable_display <- ifelse(
    !duplicated(d$block),
    d$variable,
    ""
  )
  
  
  d$n_display <- ifelse(
    is.na(d$n),
    "",
    format(
      as.integer(d$n),
      scientific = FALSE,
      trim = TRUE
    )
  )
  
  
  estimable <- (
    !d$is_reference &
      is.finite(d$estimate) &
      is.finite(d$conf.low) &
      is.finite(d$conf.high)
  )
  
  
  d$hr_display <- reference_label
  
  d$hr_display[
    !d$is_reference & !estimable
  ] <- "Not estimable"
  
  
  d$hr_display[
    estimable
  ] <- sprintf(
    paste0(
      "%.", digits,
      "f [%.", digits,
      "f, %.", digits,
      "f]"
    ),
    d$estimate[estimable],
    d$conf.low[estimable],
    d$conf.high[estimable]
  )
  
  
  d$p_display <- vapply(
    d$p.value,
    function(x) {
      
      if (is.na(x)) {
        return("")
      }
      
      .ekb_format_p(
        x,
        digits = p_digits
      )
    },
    character(1)
  )
  
  d$p_display[d$is_reference] <- ""
  
  
  # ----------------------------------------------------------
  # Vertical positions
  # ----------------------------------------------------------
  
  d$y <- 0
  
  if (nrow(d) > 1L) {
    
    for (i in 2:nrow(d)) {
      
      extra <- if (
        d$block[i] != d$block[i - 1L]
      ) {
        block_gap
      } else {
        0
      }
      
      d$y[i] <- (
        d$y[i - 1L] +
          1 +
          extra
      )
    }
  }
  
  
  d$y <- max(d$y) - d$y + 1
  
  
  # ----------------------------------------------------------
  # Forest limits and breaks
  # ----------------------------------------------------------
  
  forest_limits <- .ekb_forest_limits(
    conf.low = d$conf.low,
    conf.high = d$conf.high,
    limits = limits
  )
  
  
  forest_breaks <- .ekb_forest_breaks(
    limits = forest_limits,
    breaks = breaks
  )
  
  
  # ----------------------------------------------------------
  # Column widths
  # ----------------------------------------------------------
  
  widths <- c(
    variable = 2.4,
    level = 3.1,
    n = 0.75,
    forest = 3.2,
    hr = 2.7,
    p = 1.05
  )
  
  
  if (!is.null(column_widths)) {
    
    if (is.null(names(column_widths))) {
      stop(
        "column_widths must be a named numeric vector.",
        call. = FALSE
      )
    }
    
    known <- intersect(
      names(column_widths),
      names(widths)
    )
    
    widths[known] <- column_widths[known]
  }
  
  
  if (!isTRUE(show_n)) {
    widths["n"] <- 0
  }
  
  if (!isTRUE(show_p)) {
    widths["p"] <- 0
  }
  
  
  column_order <- c(
    "variable",
    "level",
    "n",
    "forest",
    "hr",
    "p"
  )
  
  
  widths <- widths[column_order]
  
  
  starts <- c(
    0,
    cumsum(
      widths[-length(widths)]
    )
  )
  
  names(starts) <- column_order
  
  
  total_width <- sum(widths)
  
  
  # ----------------------------------------------------------
  # Column positions
  # ----------------------------------------------------------
  
  pad <- 0.08
  
  
  x_variable <-
    starts["variable"] + pad
  
  x_level <-
    starts["level"] + pad
  
  x_n <-
    starts["n"] +
    widths["n"] / 2
  
  
  forest_left <-
    starts["forest"] + 0.15
  
  forest_right <-
    starts["forest"] +
    widths["forest"] - 0.15
  
  
  # Right-align numerical columns
  x_hr <-
    starts["hr"] +
    widths["hr"] - pad
  
  x_p <-
    starts["p"] +
    widths["p"] - pad
  
  
  # ----------------------------------------------------------
  # Transform HR onto forest-panel coordinates
  # ----------------------------------------------------------
  
  map_hr <- function(x) {
    
    forest_left +
      (
        log(x) -
          log(forest_limits[1])
      ) /
      (
        log(forest_limits[2]) -
          log(forest_limits[1])
      ) *
      (
        forest_right -
          forest_left
      )
  }
  
  
  reference_x <- map_hr(reference)
  
  
  # ----------------------------------------------------------
  # Plot coordinates
  # ----------------------------------------------------------
  
  header_y <-
    max(d$y) + 0.85
  
  header_line_y <-
    max(d$y) + 0.42
  
  bottom_line_y <-
    min(d$y) - 0.48
  
  axis_y <-
    min(d$y) - 0.75
  
  tick_length <- 0.16
  
  tick_y <-
    axis_y - 0.38
  
  use_favors <- (
    isTRUE(show_favors) &&
      !is.null(favors_left) &&
      !is.null(favors_right)
  )
  
  if (use_favors) {
    favors_arrow_y <-
      tick_y - 0.46
    
    favors_text_y <-
      favors_arrow_y - 0.48
    
    xlab_y <-
      favors_text_y - 0.66
  } else {
    xlab_y <-
      tick_y - 0.55
  }
  
  
  
  use_title <-
    !is.null(title) &&
    nzchar(as.character(title))
  
  title_y <- if (use_title) {
    header_y + 0.72
  } else {
    header_y
  }
  
  # ----------------------------------------------------------
  # Forest data
  # ----------------------------------------------------------
  
  points <- d[
    estimable,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(points)) {
    
    points$estimate_x <- map_hr(
      pmin(
        pmax(
          points$estimate,
          forest_limits[1]
        ),
        forest_limits[2]
      )
    )
    
    points$low_x <- map_hr(
      pmin(
        pmax(
          points$conf.low,
          forest_limits[1]
        ),
        forest_limits[2]
      )
    )
    
    points$high_x <- map_hr(
      pmin(
        pmax(
          points$conf.high,
          forest_limits[1]
        ),
        forest_limits[2]
      )
    )
  }
  
  
  tick_df <- data.frame(
    value = forest_breaks,
    x = map_hr(forest_breaks),
    stringsAsFactors = FALSE
  )
  
  
  tick_df$label <- format(
    tick_df$value,
    scientific = FALSE,
    trim = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Base plot
  # ----------------------------------------------------------
  
  p <- ggplot2::ggplot()
  
  
  # ----------------------------------------------------------
  # Alternating row shading
  # ----------------------------------------------------------
  
  if (isTRUE(shade)) {
    
    shade_rows <- seq(
      1L,
      nrow(d),
      by = 2L
    )
    
    
    shade_df <- d[
      shade_rows,
      ,
      drop = FALSE
    ]
    
    
    p <- p +
      ggplot2::geom_rect(
        data = shade_df,
        ggplot2::aes(
          xmin = 0,
          xmax = total_width,
          ymin = .data$y - 0.43,
          ymax = .data$y + 0.43
        ),
        inherit.aes = FALSE,
        fill = "grey94",
        colour = NA
      )
  }
  
  
  # ----------------------------------------------------------
  # Header and bottom rules
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = total_width,
        y = header_line_y,
        yend = header_line_y
      ),
      linewidth = 0.5,
      colour = "black"
    ) +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = total_width,
        y = bottom_line_y,
        yend = bottom_line_y
      ),
      linewidth = 0.45,
      colour = "black"
    )
  
  
  # ----------------------------------------------------------
  # Left table columns
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(
        y = .data$y,
        label = .data$variable_display
      ),
      x = x_variable,
      hjust = 0,
      fontface = "bold",
      size = text_size
    ) +
    
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(
        y = .data$y,
        label = .data$level
      ),
      x = x_level,
      hjust = 0,
      size = text_size
    )
  
  
  if (isTRUE(show_n)) {
    
    p <- p +
      ggplot2::geom_text(
        data = d,
        ggplot2::aes(
          y = .data$y,
          label = .data$n_display
        ),
        x = x_n,
        hjust = 0.5,
        size = text_size
      )
  }
  
  
  # ----------------------------------------------------------
  # Reference line
  # ----------------------------------------------------------
  
  p <- p +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = reference_x,
        xend = reference_x,
        y = bottom_line_y,
        yend = header_line_y
      ),
      linewidth = reference_linewidth,
      linetype = "dotted",
      colour = "grey25"
    )
  
  
  # ----------------------------------------------------------
  # Confidence intervals
  # ----------------------------------------------------------
  
  if (nrow(points)) {
    
    cap_height <- 0.13
    
    
    p <- p +
      
      ggplot2::geom_segment(
        data = points,
        ggplot2::aes(
          x = .data$low_x,
          xend = .data$high_x,
          y = .data$y,
          yend = .data$y
        ),
        linewidth = ci_linewidth,
        colour = "black"
      ) +
      
      # left CI cap
      ggplot2::geom_segment(
        data = points,
        ggplot2::aes(
          x = .data$low_x,
          xend = .data$low_x,
          y = .data$y - cap_height,
          yend = .data$y + cap_height
        ),
        linewidth = ci_linewidth,
        colour = "black"
      ) +
      
      # right CI cap
      ggplot2::geom_segment(
        data = points,
        ggplot2::aes(
          x = .data$high_x,
          xend = .data$high_x,
          y = .data$y - cap_height,
          yend = .data$y + cap_height
        ),
        linewidth = ci_linewidth,
        colour = "black"
      ) +
      
      # square estimate
      ggplot2::geom_point(
        data = points,
        ggplot2::aes(
          x = .data$estimate_x,
          y = .data$y
        ),
        size = point_size,
        shape = 15,
        colour = "black"
      )
  }
  
  
  # ----------------------------------------------------------
  # Right table columns
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(
        y = .data$y,
        label = .data$hr_display
      ),
      x = x_hr,
      hjust = 1,
      size = text_size
    )
  
  
  if (isTRUE(show_p)) {
    
    p <- p +
      ggplot2::geom_text(
        data = d,
        ggplot2::aes(
          y = .data$y,
          label = .data$p_display
        ),
        x = x_p,
        hjust = 1,
        size = text_size
      )
  }
  
  
  # ----------------------------------------------------------
  # Headers
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::annotate(
      "text",
      x = x_variable,
      y = header_y,
      label = variable_heading,
      hjust = 0,
      fontface = "bold",
      size = header_size
    ) +
    
    ggplot2::annotate(
      "text",
      x = x_level,
      y = header_y,
      label = level_heading,
      hjust = 0,
      fontface = "bold",
      size = header_size
    )
  
  
  if (isTRUE(show_n)) {
    
    p <- p +
      ggplot2::annotate(
        "text",
        x = x_n,
        y = header_y,
        label = n_heading,
        hjust = 0.5,
        fontface = "bold",
        size = header_size
      )
  }
  
  
  p <- p +
    ggplot2::annotate(
      "text",
      x = x_hr,
      y = header_y,
      label = hr_heading,
      hjust = 1,
      fontface = "bold",
      size = header_size
    )
  
  
  if (isTRUE(show_p)) {
    
    p <- p +
      ggplot2::annotate(
        "text",
        x = x_p,
        y = header_y,
        label = p_heading,
        hjust = 1,
        fontface = "bold",
        size = header_size
      )
  }
  
  
  # ----------------------------------------------------------
  # Forest x-axis
  # ----------------------------------------------------------
  
  p <- p +
    
    ggplot2::geom_segment(
      ggplot2::aes(
        x = forest_left,
        xend = forest_right,
        y = axis_y,
        yend = axis_y
      ),
      linewidth = 0.45,
      colour = "black"
    ) +
    
    ggplot2::geom_segment(
      data = tick_df,
      ggplot2::aes(
        x = .data$x,
        xend = .data$x,
        y = axis_y,
        yend = axis_y - tick_length
      ),
      linewidth = 0.45,
      colour = "black"
    ) +
    
    ggplot2::geom_text(
      data = tick_df,
      ggplot2::aes(
        x = .data$x,
        y = tick_y,
        label = .data$label
      ),
      size = text_size,
      hjust = 0.5
    )
  
  
  # ----------------------------------------------------------
  # Plot title
  # ----------------------------------------------------------
  
  if (use_title) {
    
    p <- p +
      ggplot2::annotate(
        "text",
        x = total_width / 2,
        y = title_y,
        label = title,
        hjust = 0.5,
        fontface = "bold",
        size = header_size * 1.12
      )
  }
  
  
  # ----------------------------------------------------------
  # Treatment direction / favors annotations
  # ----------------------------------------------------------
  
  if (use_favors) {
    
    forest_span <-
      forest_right - forest_left
    
    reference_gap <-
      0.04 * forest_span
    
    outer_gap <-
      0.03 * forest_span
    
    # HR < 1: lower hazard for comparator treatment
    p <- p +
      ggplot2::annotate(
        "segment",
        x = reference_x - reference_gap,
        xend = forest_left + outer_gap,
        y = favors_arrow_y,
        yend = favors_arrow_y,
        linewidth = 0.4,
        colour = "black",
        arrow = grid::arrow(
          length = grid::unit(0.07, "in"),
          type = "closed"
        )
      ) +
      ggplot2::annotate(
        "text",
        x = (reference_x + forest_left) / 2,
        y = favors_text_y,
        label = paste0(favors_left, " better"),
        hjust = 0.5,
        size = text_size
      )
    
    # HR > 1: lower hazard for reference treatment
    p <- p +
      ggplot2::annotate(
        "segment",
        x = reference_x + reference_gap,
        xend = forest_right - outer_gap,
        y = favors_arrow_y,
        yend = favors_arrow_y,
        linewidth = 0.4,
        colour = "black",
        arrow = grid::arrow(
          length = grid::unit(0.07, "in"),
          type = "closed"
        )
      ) +
      ggplot2::annotate(
        "text",
        x = (reference_x + forest_right) / 2,
        y = favors_text_y,
        label = paste0(favors_right, " better"),
        hjust = 0.5,
        size = text_size
      )
  }
  
  
  if (!is.null(xlab)) {
    
    p <- p +
      ggplot2::annotate(
        "text",
        x = (
          forest_left +
            forest_right
        ) / 2,
        y = xlab_y,
        label = xlab,
        hjust = 0.5,
        size = axis_title_size
      )
  }
  
  
  # ----------------------------------------------------------
  # Final theme
  # ----------------------------------------------------------
  
  p +
    ggplot2::coord_cartesian(
      xlim = c(
        0,
        total_width
      ),
      ylim = c(
        xlab_y - 0.35,
        if (use_title) title_y + 0.38 else header_y + 0.35
      ),
      clip = "off"
    ) +
    ggplot2::theme_void(
      base_size = base_size
    ) +
    ggplot2::theme(
      plot.background =
        ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
      
      panel.background =
        ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
      
      plot.margin =
        ggplot2::margin(
          t = 8,
          r = 8,
          b = 8,
          l = 8
        )
    )
}


# ============================================================
# Prepare Cox model for publication-style forest plot
# ============================================================

.ekb_cox_forest_data <- function(
    object,
    labels = NULL,
    level_labels = NULL,
    show_reference = TRUE,
    continuous_label = "Per unit increase",
    conf.level = 0.95
) {
  
  if (!inherits(object, "ekb_cox_fit")) {
    stop(
      "object must be created by fit_cox().",
      call. = FALSE
    )
  }
  
  
  z <- as.data.frame(
    tidy_cox(
      object,
      conf.level = conf.level
    )
  )
  
  
  fit <- object$fit
  covariates <- object$covariates
  xlevels <- fit$xlevels
  model_data <- fit$model
  
  out <- list()
  k <- 1L
  cursor <- 1L
  
  
  for (variable in covariates) {
    
    variable_label <- .ekb_forest_label(
      variable,
      labels
    )
    
    
    levels_variable <- xlevels[[variable]]
    
    
    # --------------------------------------------------------
    # Categorical variable
    # --------------------------------------------------------
    
    if (
      !is.null(levels_variable) &&
      length(levels_variable) >= 2L
    ) {
      
      n_coef <-
        length(levels_variable) - 1L
      
      
      idx <- seq.int(
        cursor,
        length.out = n_coef
      )
      
      
      if (max(idx) > nrow(z)) {
        stop(
          sprintf(
            "Could not reconstruct coefficient rows for '%s'.",
            variable
          ),
          call. = FALSE
        )
      }
      
      
      zz <- z[
        idx,
        ,
        drop = FALSE
      ]
      
      
      cursor <-
        cursor + n_coef
      
      
      block_rows <- list()
      j <- 1L
      
      
      if (isTRUE(show_reference)) {
        
        reference_level <-
          levels_variable[1]
        
        
        block_rows[[j]] <- data.frame(
          block = variable,
          variable = variable_label,
          level = .ekb_forest_level_label(
            variable,
            reference_level,
            level_labels
          ),
          estimate = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          n = sum(
            as.character(model_data[[variable]]) ==
              as.character(reference_level),
            na.rm = TRUE
          ),
          is_reference = TRUE,
          stringsAsFactors = FALSE
        )
        
        j <- j + 1L
      }
      
      
      for (i in seq_len(n_coef)) {
        
        level_i <-
          levels_variable[i + 1L]
        
        
        block_rows[[j]] <- data.frame(
          block = variable,
          variable = variable_label,
          level = .ekb_forest_level_label(
            variable,
            level_i,
            level_labels
          ),
          estimate = zz$hr[i],
          conf.low = zz$hr_conf.low[i],
          conf.high = zz$hr_conf.high[i],
          p.value = zz$p.value[i],
          n = sum(
            as.character(model_data[[variable]]) ==
              as.character(level_i),
            na.rm = TRUE
          ),
          is_reference = FALSE,
          stringsAsFactors = FALSE
        )
        
        j <- j + 1L
      }
      
      
      out[[k]] <-
        do.call(
          rbind,
          block_rows
        )
      
      k <- k + 1L
      
      next
    }
    
    
    # --------------------------------------------------------
    # Continuous variable
    # --------------------------------------------------------
    
    if (cursor > nrow(z)) {
      stop(
        sprintf(
          "Could not reconstruct coefficient row for '%s'.",
          variable
        ),
        call. = FALSE
      )
    }
    
    
    zz <- z[
      cursor,
      ,
      drop = FALSE
    ]
    
    
    cursor <-
      cursor + 1L
    
    
    out[[k]] <- data.frame(
      block = variable,
      variable = variable_label,
      level = continuous_label,
      estimate = zz$hr,
      conf.low = zz$hr_conf.low,
      conf.high = zz$hr_conf.high,
      p.value = zz$p.value,
      n = fit$n,
      is_reference = FALSE,
      stringsAsFactors = FALSE
    )
    
    
    k <- k + 1L
  }
  
  
  do.call(
    rbind,
    out
  )
}


# ============================================================
# Prepare multiply-imputed Cox model for publication-style forest plot
# ============================================================

.ekb_mi_mean_n <- function(
    object,
    variable,
    level = NULL
) {
  
  .ekb_require("mice")
  
  if (is.null(object$imputation)) {
    stop(
      "MI object does not contain an imputation component.",
      call. = FALSE
    )
  }
  
  mi <- if (inherits(object$imputation, "ekb_imputation")) {
    object$imputation$mids
  } else {
    object$imputation
  }
  
  if (!inherits(mi, "mids")) {
    stop(
      "MI object does not contain a valid mice::mids object.",
      call. = FALSE
    )
  }
  
  counts <- vapply(
    seq_len(mi$m),
    function(i) {
      
      d <- mice::complete(
        mi,
        i
      )
      
      if (is.null(level)) {
        return(
          sum(!is.na(d[[variable]]))
        )
      }
      
      sum(
        as.character(d[[variable]]) ==
          as.character(level),
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  # N in MI forests is the rounded mean category size across
  # completed datasets. It is descriptive, not Rubin-pooled.
  as.integer(
    round(
      mean(counts)
    )
  )
}


.ekb_mi_cox_covariates <- function(object) {
  
  if (inherits(object, "ekb_mi_cox")) {
    return(object$covariates)
  }
  
  if (inherits(object, "ekb_mi_iptw_cox")) {
    return(
      unique(
        c(
          object$treatment,
          setdiff(
            object$outcome_covariates %||% character(),
            object$treatment
          )
        )
      )
    )
  }
  
  stop(
    paste(
      "object must be created by fit_mi_cox()."
    ),
    call. = FALSE
  )
}


.ekb_mi_cox_forest_data <- function(
    object,
    labels = NULL,
    level_labels = NULL,
    show_reference = TRUE,
    continuous_label = "Per unit increase"
) {
  
  if (!inherits(
    object,
    c(
      "ekb_mi_cox",
      "ekb_mi_iptw_cox"
    )
  )) {
    stop(
      "object must be created by fit_mi_cox().",
      call. = FALSE
    )
  }
  
  if (
    is.null(object$pooled_summary) ||
    !nrow(as.data.frame(object$pooled_summary))
  ) {
    stop(
      "MI object contains no pooled Cox results.",
      call. = FALSE
    )
  }
  
  if (
    is.null(object$fits) ||
    !length(object$fits)
  ) {
    stop(
      "MI object contains no fitted Cox models.",
      call. = FALSE
    )
  }
  
  z <- as.data.frame(
    object$pooled_summary
  )
  
  required <- c(
    "term",
    "estimate",
    "conf.low",
    "conf.high",
    "p.value"
  )
  
  .ekb_assert_columns(
    z,
    required
  )
  
  fit <- object$fits[[1L]]
  covariates <- .ekb_mi_cox_covariates(object)
  xlevels <- fit$xlevels
  
  out <- list()
  k <- 1L
  cursor <- 1L
  
  for (variable in covariates) {
    
    variable_label <- .ekb_forest_label(
      variable,
      labels
    )
    
    levels_variable <- xlevels[[variable]]
    
    # --------------------------------------------------------
    # Categorical variable
    # --------------------------------------------------------
    
    if (
      !is.null(levels_variable) &&
      length(levels_variable) >= 2L
    ) {
      
      n_coef <-
        length(levels_variable) - 1L
      
      idx <- seq.int(
        cursor,
        length.out = n_coef
      )
      
      if (max(idx) > nrow(z)) {
        stop(
          sprintf(
            "Could not reconstruct pooled coefficient rows for '%s'.",
            variable
          ),
          call. = FALSE
        )
      }
      
      zz <- z[
        idx,
        ,
        drop = FALSE
      ]
      
      cursor <-
        cursor + n_coef
      
      block_rows <- list()
      j <- 1L
      
      if (isTRUE(show_reference)) {
        
        reference_level <-
          levels_variable[1L]
        
        block_rows[[j]] <- data.frame(
          block = variable,
          variable = variable_label,
          level = .ekb_forest_level_label(
            variable,
            reference_level,
            level_labels
          ),
          estimate = NA_real_,
          conf.low = NA_real_,
          conf.high = NA_real_,
          p.value = NA_real_,
          n = .ekb_mi_mean_n(
            object,
            variable,
            reference_level
          ),
          is_reference = TRUE,
          stringsAsFactors = FALSE
        )
        
        j <- j + 1L
      }
      
      for (i in seq_len(n_coef)) {
        
        level_i <-
          levels_variable[i + 1L]
        
        block_rows[[j]] <- data.frame(
          block = variable,
          variable = variable_label,
          level = .ekb_forest_level_label(
            variable,
            level_i,
            level_labels
          ),
          estimate = zz$estimate[i],
          conf.low = zz$conf.low[i],
          conf.high = zz$conf.high[i],
          p.value = zz$p.value[i],
          n = .ekb_mi_mean_n(
            object,
            variable,
            level_i
          ),
          is_reference = FALSE,
          stringsAsFactors = FALSE
        )
        
        j <- j + 1L
      }
      
      out[[k]] <-
        do.call(
          rbind,
          block_rows
        )
      
      k <- k + 1L
      
      next
    }
    
    # --------------------------------------------------------
    # Continuous variable
    # --------------------------------------------------------
    
    if (cursor > nrow(z)) {
      stop(
        sprintf(
          "Could not reconstruct pooled coefficient row for '%s'.",
          variable
        ),
        call. = FALSE
      )
    }
    
    zz <- z[
      cursor,
      ,
      drop = FALSE
    ]
    
    cursor <-
      cursor + 1L
    
    out[[k]] <- data.frame(
      block = variable,
      variable = variable_label,
      level = continuous_label,
      estimate = zz$estimate,
      conf.low = zz$conf.low,
      conf.high = zz$conf.high,
      p.value = zz$p.value,
      n = .ekb_mi_mean_n(
        object,
        variable
      ),
      is_reference = FALSE,
      stringsAsFactors = FALSE
    )
    
    k <- k + 1L
  }
  
  do.call(
    rbind,
    out
  )
}


# ============================================================
# Publication-style Cox forest plot
# ============================================================

#' Plot Cox model effects
#'
#' Plots ekb_cox_fit or ekb_mi_cox effects with categorical reference rows.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param level_labels Optional named list of named character vectors mapping variable levels to display labels.
#' @param show_reference Whether to show reference-category rows.
#' @param show_n Whether to show observation counts.
#' @param show_p Whether to display p-values.
#' @param title Optional plot title.
#' @param xlab Horizontal axis label.
#' @param limits Optional positive pair of hazard-ratio axis limits.
#' @param breaks Optional positive hazard-ratio tick positions.
#' @param digits Number of decimal places for estimates.
#' @param p_digits Number of decimal places for p-values.
#' @param base_size Base plot font size in points.
#' @param point_size Point size in the forest panel.
#' @param block_gap Vertical spacing between variable blocks.
#' @param column_widths Optional relative widths of the forest-plot columns.
#' @param conf.level Confidence level between zero and one.
#' @return A ggplot forest plot.
#' @details Reference rows and counts are derived from the fitted model data. MI counts are averaged across completed datasets. The pooled MI intervals stored in the model are used; conf.level controls extraction for ordinary Cox fits.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit <- fit_cox(d, "time", "event", c("treatment", "age"))
#' if (requireNamespace("patchwork", quietly = TRUE)) plot_cox_forest(fit)
plot_cox_forest <- function(
    object,
    labels = NULL,
    level_labels = NULL,
    show_reference = TRUE,
    show_n = TRUE,
    show_p = TRUE,
    title = NULL,
    xlab = "Hazard ratio",
    limits = NULL,
    breaks = NULL,
    digits = 2,
    p_digits = 3,
    base_size = 10,
    point_size = 2.2,
    block_gap = 0.35,
    column_widths = NULL,
    conf.level = 0.95
) {
  
  if (inherits(object, "ekb_cox_fit")) {
    
    d <- .ekb_cox_forest_data(
      object = object,
      labels = labels,
      level_labels = level_labels,
      show_reference = show_reference,
      conf.level = conf.level
    )
    
  } else if (inherits(
    object,
    c(
      "ekb_mi_cox",
      "ekb_mi_iptw_cox"
    )
  )) {
    
    d <- .ekb_mi_cox_forest_data(
      object = object,
      labels = labels,
      level_labels = level_labels,
      show_reference = show_reference
    )
    
  } else {
    
    stop(
      paste(
        "object must be created by fit_cox() or fit_mi_cox()."
      ),
      call. = FALSE
    )
  }
  
  .ekb_plot_forest(
    data = d,
    show_n = show_n,
    show_p = show_p,
    variable_heading = "Variable",
    level_heading = "Category",
    n_heading = "N",
    hr_heading = "HR (95% CI)",
    p_heading = "P value",
    title = title,
    xlab = xlab,
    limits = limits,
    breaks = breaks,
    digits = digits,
    p_digits = p_digits,
    base_size = base_size,
    point_size = point_size,
    block_gap = block_gap,
    column_widths = column_widths
  )
}


# ============================================================
# Publication-style subgroup forest plot
# ============================================================

#' Plot subgroup treatment effects
#'
#' Plots ekb_subgroup_cox or ekb_subgroup_mi_cox results.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param labels Optional named character vector mapping variable names to display labels.
#' @param level_labels Optional named list of named character vectors mapping variable levels to display labels.
#' @param show_n Whether to show observation counts.
#' @param p_type Display within-level treatment-effect p-values (default), available interaction p-values, or none.
#' @param show_favors Whether to label the direction of the treatment effect.
#' @param treatment_labels Optional named character vector mapping treatment levels to display labels.
#' @param lower_hr_better Whether lower hazard is interpreted as favourable for the comparator.
#' @param title Optional plot title.
#' @param xlab Horizontal axis label.
#' @param limits Optional positive pair of hazard-ratio axis limits.
#' @param breaks Optional positive hazard-ratio tick positions.
#' @param digits Number of decimal places for estimates.
#' @param p_digits Number of decimal places for p-values.
#' @param base_size Base plot font size in points.
#' @param point_size Point size in the forest panel.
#' @param block_gap Vertical spacing between variable blocks.
#' @param column_widths Optional relative widths of the forest-plot columns.
#' @return A ggplot forest plot.
#' @details The default p_type = "effect" displays treatment-effect p-values within each level. Unavailable interaction tests warn and fall back to effect p-values. N is an observation count, not weighted N; MI counts are averaged over completed datasets.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' sg <- subgroup_cox(d, "time", "event", "treatment", "sex")
#' if (requireNamespace("patchwork", quietly = TRUE)) plot_subgroup_forest(sg)
plot_subgroup_forest <- function(
    object,
    labels = NULL,
    level_labels = NULL,
    show_n = TRUE,
    p_type = c(
      "effect",
      "interaction",
      "none"
    ),
    show_favors = TRUE,
    treatment_labels = NULL,
    lower_hr_better = TRUE,
    title = NULL,
    xlab = "Hazard ratio",
    limits = NULL,
    breaks = NULL,
    digits = 2,
    p_digits = 3,
    base_size = 10,
    point_size = 2.2,
    block_gap = 0.35,
    column_widths = NULL
) {
  
  if (!inherits(
    object,
    c(
      "ekb_subgroup_cox",
      "ekb_subgroup_mi_cox"
    )
  )) {
    stop(
      paste(
        "object must be created by",
        "subgroup_cox() or subgroup_mi_cox()."
      ),
      call. = FALSE
    )
  }
  
  r <- as.data.frame(
    object$results
  )
  
  if (!nrow(r)) {
    stop(
      "No subgroup results to plot.",
      call. = FALSE
    )
  }
  
  p_type <- match.arg(
    p_type
  )
  
  # ----------------------------------------------------------
  # Interaction availability
  # ----------------------------------------------------------
  
  interaction_available <- (
    !is.null(object$interaction_tests) &&
      nrow(
        as.data.frame(
          object$interaction_tests
        )
      ) > 0L
  )
  
  if (
    p_type == "interaction" &&
    !interaction_available
  ) {
    warning(
      paste(
        "No interaction tests are available.",
        "Falling back to subgroup-specific p-values."
      ),
      call. = FALSE
    )
    p_type <- "effect"
  }
  
  # ----------------------------------------------------------
  # Preserve original subgroup order
  # ----------------------------------------------------------
  
  subgroup_order <- unique(
    as.character(
      r$subgroup[
        r$subgroup != "Overall"
      ]
    )
  )
  
  overall <- r[
    r$subgroup == "Overall",
    ,
    drop = FALSE
  ]
  
  blocks <- list()
  k <- 1L
  
  # ----------------------------------------------------------
  # Overall effect first
  # ----------------------------------------------------------
  
  if (nrow(overall)) {
    blocks[[k]] <- data.frame(
      block = "Overall",
      variable = "Overall",
      level = "",
      estimate = overall$hr[1],
      conf.low = overall$conf.low[1],
      conf.high = overall$conf.high[1],
      p.value = if (
        p_type == "effect"
      ) {
        overall$p.value[1]
      } else {
        NA_real_
      },
      n = overall$n[1],
      is_reference = FALSE,
      stringsAsFactors = FALSE
    )
    
    k <- k + 1L
  }
  
  # ----------------------------------------------------------
  # Interaction p-values
  # ----------------------------------------------------------
  
  interaction_table <- NULL
  
  if (interaction_available) {
    interaction_table <- as.data.frame(
      object$interaction_tests
    )
  }
  
  # ----------------------------------------------------------
  # Subgroup blocks
  # ----------------------------------------------------------
  
  for (subgroup in subgroup_order) {
    
    rr <- r[
      r$subgroup == subgroup,
      ,
      drop = FALSE
    ]
    
    subgroup_label <- .ekb_forest_label(
      subgroup,
      labels
    )
    
    p_interaction <- NA_real_
    
    if (
      p_type == "interaction" &&
      !is.null(interaction_table)
    ) {
      hit <- interaction_table[
        interaction_table$subgroup == subgroup,
        ,
        drop = FALSE
      ]
      
      if (nrow(hit)) {
        p_interaction <-
          hit$interaction_p[1]
      }
    }
    
    for (i in seq_len(nrow(rr))) {
      
      level_i <- as.character(
        rr$level[i]
      )
      
      blocks[[k]] <- data.frame(
        block = subgroup,
        variable = subgroup_label,
        level = .ekb_forest_level_label(
          subgroup,
          level_i,
          level_labels
        ),
        estimate = rr$hr[i],
        conf.low = rr$conf.low[i],
        conf.high = rr$conf.high[i],
        p.value = if (
          p_type == "interaction"
        ) {
          if (i == 1L) {
            p_interaction
          } else {
            NA_real_
          }
        } else if (
          p_type == "effect"
        ) {
          rr$p.value[i]
        } else {
          NA_real_
        },
        n = rr$n[i],
        is_reference = FALSE,
        stringsAsFactors = FALSE
      )
      
      k <- k + 1L
    }
  }
  
  d <- do.call(
    rbind,
    blocks
  )
  
  # ----------------------------------------------------------
  # P-value heading
  # ----------------------------------------------------------
  
  show_p <-
    p_type != "none"
  
  p_heading <- if (
    p_type == "interaction"
  ) {
    "P interaction"
  } else {
    "P value"
  }
  
  # ----------------------------------------------------------
  # Treatment labels for direction arrows
  # ----------------------------------------------------------
  
  treatment_reference <-
    object$treatment_reference
  
  treatment_comparator <-
    object$treatment_comparator
  
  if (
    is.null(treatment_reference) ||
    is.null(treatment_comparator)
  ) {
    show_favors <- FALSE
  }
  
  display_reference <-
    treatment_reference
  
  display_comparator <-
    treatment_comparator
  
  if (!is.null(treatment_labels)) {
    
    if (is.null(names(treatment_labels))) {
      stop(
        "treatment_labels must be a named character vector.",
        call. = FALSE
      )
    }
    
    if (
      !is.null(display_reference) &&
      display_reference %in%
      names(treatment_labels)
    ) {
      display_reference <-
        treatment_labels[[display_reference]]
    }
    
    if (
      !is.null(display_comparator) &&
      display_comparator %in%
      names(treatment_labels)
    ) {
      display_comparator <-
        treatment_labels[[display_comparator]]
    }
  }
  
  # Cox coefficient is comparator vs reference.
  # Therefore HR < 1 favours the comparator when a lower
  # event hazard represents a better outcome.
  if (isTRUE(lower_hr_better)) {
    favors_left <-
      display_comparator
    favors_right <-
      display_reference
  } else {
    favors_left <-
      display_reference
    favors_right <-
      display_comparator
  }
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  .ekb_plot_forest(
    data = d,
    show_n = show_n,
    show_p = show_p,
    variable_heading = "Subgroup",
    level_heading = "Category",
    n_heading = "N",
    hr_heading = "HR (95% CI)",
    p_heading = p_heading,
    title = title,
    xlab = xlab,
    limits = limits,
    breaks = breaks,
    digits = digits,
    p_digits = p_digits,
    base_size = base_size,
    point_size = point_size,
    block_gap = block_gap,
    column_widths = column_widths,
    show_favors = show_favors,
    favors_left = favors_left,
    favors_right = favors_right
  )
}




