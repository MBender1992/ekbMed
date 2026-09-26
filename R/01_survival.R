# Survival / Kaplan-Meier functions

#' Fit Kaplan-Meier survival curves
#'
#' Fits unweighted or case-weighted curves through survival::survfit.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param group Optional name of the grouping column.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @param conf.type Confidence-interval transformation passed to [survival::survfit()].
#' @param na.action Missing-data handler passed to [survival::survfit()].
#' @param conf.level Confidence level between zero and one.
#' @return An `ekb_survival_fit` list containing fit, input data, endpoint names and resolved numeric weights.
#' @details Weights are reused without estimating a propensity model. The survfit backend controls weighted variance and confidence intervals. Summary counts and plot risk-table counts are raw observations; weighted curve estimates are distinct from these counts.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' fit_survival(d, "time", "event", "treatment")
fit_survival <- function(data,
                         time,
                         event,
                         group = NULL,
                         weights = NULL,
                         conf.type = "log-log",
                         na.action = stats::na.omit,
                         conf.level = 0.95) {
  .ekb_require("survival")
  .ekb_assert_columns(data, c(time, event, group))

  dat <- as.data.frame(data)
  dat[[event]] <- .ekb_binary_event(dat[[event]], event)
  if (!is.null(group)) dat[[group]] <- droplevels(factor(dat[[group]]))
  .ekb_validate_time(dat, time)
  w <- .ekb_resolve_weights(dat, weights)

  formula <- .ekb_surv_formula(time, event, group)
  args <- list(
    formula = formula,
    data = dat,
    conf.type = conf.type,
    conf.int = conf.level,
    na.action = na.action,
    model = TRUE
  )
  if (!is.null(w)) args$weights <- w
  fit <- do.call(survival::survfit, args)

  structure(
    list(
      fit = fit,
      data = dat,
      time = time,
      event = event,
      group = group,
      weights = w,
      conf.type = conf.type,
      na.action = na.action,
      call = match.call()
    ),
    class = "ekb_survival_fit"
  )
}

#' Compare survival curves
#'
#' Runs a standard log-rank test or the RISCA weighted log-rank test.
#'
#' @param data A data frame. Variable arguments are column names supplied as character strings.
#' @param time Name of the non-negative follow-up or stop-time column. Use months for month-labelled tables.
#' @param event Name of the event column, coded 0/1 or logical; 1 means the event occurred.
#' @param group Optional name of the grouping column.
#' @param weights NULL, a positive numeric vector aligned to the input rows, a weight-column name, or an `ekb_weighting` object. Do not reorder or filter rows after estimating weights without realigning them.
#' @return A list with statistic, df, p.value and method.
#' @details Without weights uses survival::survdiff. With weights uses RISCA::ipw.log.rank and requires two observed groups. Rows missing endpoint, group or weights are removed explicitly.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' survival_logrank(d, "time", "event", "treatment")
survival_logrank <- function(data,
                             time,
                             event,
                             group,
                             weights = NULL) {
  .ekb_require("survival")
  .ekb_assert_columns(data, c(time, event, group))
  dat <- as.data.frame(data)
  dat[[event]] <- .ekb_binary_event(dat[[event]], event)
  dat[[group]] <- droplevels(factor(dat[[group]]))
  .ekb_validate_time(dat, time)
  w <- .ekb_resolve_weights(dat, weights)

  keep <- stats::complete.cases(dat[, c(time, event, group), drop = FALSE])
  if (!is.null(w)) keep <- keep & !is.na(w)
  dat <- dat[keep, , drop = FALSE]
  if (!is.null(w)) w <- w[keep]

  dat[[group]] <- droplevels(dat[[group]])
  if (nlevels(dat[[group]]) < 2L) stop("At least two groups are required.", call. = FALSE)

  if (is.null(w)) {
    formula <- .ekb_surv_formula(time, event, group)
    test <- survival::survdiff(formula, data = dat)
    df <- length(test$n) - 1L
    p <- stats::pchisq(test$chisq, df = df, lower.tail = FALSE)
    return(list(statistic = unname(test$chisq), df = df, p.value = p, method = "Log-rank test"))
  }

  if (nlevels(dat[[group]]) != 2L) {
    stop("Weighted log-rank testing via RISCA supports exactly two groups.", call. = FALSE)
  }
  .ekb_require("RISCA")
  g <- as.integer(dat[[group]] == levels(dat[[group]])[2L])
  test <- RISCA::ipw.log.rank(
    times = dat[[time]],
    failures = dat[[event]],
    variable = g,
    weights = w
  )
  list(
    statistic = unname(test$statistic),
    df = 1L,
    p.value = unname(test$p.value),
    method = "IPTW-adjusted log-rank test (RISCA::ipw.log.rank)"
  )
}

.ekb_survival_raw_counts <- function(object) {
  dat <- object$data
  keep <- rep(TRUE, nrow(dat))
  if (!is.null(object$fit$na.action)) keep[as.integer(object$fit$na.action)] <- FALSE
  if (is.null(object$group)) {
    return(data.frame(
      group = "Overall",
      n = sum(keep),
      events = sum(dat[[object$event]][keep] == 1L),
      stringsAsFactors = FALSE
    ))
  }
  dat <- dat[keep, , drop = FALSE]
  lev <- levels(droplevels(dat[[object$group]]))
  do.call(rbind, lapply(lev, function(z) {
    idx <- dat[[object$group]] == z
    data.frame(
      group = z,
      n = sum(idx),
      events = sum(dat[[object$event]][idx] == 1L),
      stringsAsFactors = FALSE
    )
  }))
}

#' Summarize survival medians and landmarks
#'
#' Summarizes an ekb_survival_fit object without extrapolation by default.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param landmarks Numeric follow-up times in the same unit as the time column. NULL or `numeric()` requests no landmarks.
#' @param conf.level Confidence level between zero and one.
#' @param extend Whether to extend estimates beyond observed follow-up; FALSE avoids extrapolation.
#' @return A list containing median, landmarks, raw_counts, weighted, conf.level and extend.
#' @details Median n and events are raw counts, including for weighted fits. Landmark n.risk values originate from survfit and can be weighted. A different conf.level refits the same survival specification. Times beyond follow-up are omitted unless extend is TRUE.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' s <- fit_survival(d, "time", "event", "treatment")
#' summarize_survival(s, landmarks = c(12, 24))
summarize_survival <- function(object,
                               landmarks = c(12, 24, 36, 48),
                               conf.level = 0.95,
                               extend = FALSE) {
  if (!inherits(object, "ekb_survival_fit")) {
    stop("object must be created by fit_survival().", call. = FALSE)
  }
  fit <- object$fit
  if (!isTRUE(all.equal(fit$conf.int, conf.level))) {
    fit <- fit_survival(object$data, object$time, object$event, object$group,
                         weights = object$weights, conf.type = object$conf.type,
                         na.action = object$na.action %||% stats::na.omit,
                         conf.level = conf.level)$fit
  }
  tab <- summary(fit)$table

  if (is.null(dim(tab))) {
    tab <- matrix(tab, nrow = 1L, dimnames = list("Overall", names(tab)))
  }
  tab <- as.data.frame(tab, check.names = FALSE)
  tab$group <- rownames(tab)
  rownames(tab) <- NULL
  tab$group <- sub("^[^=]+=", "", tab$group)

  med_col <- if ("median" %in% names(tab)) "median" else grep("median", names(tab), value = TRUE)[1]
  low_col <- grep("LCL", names(tab), value = TRUE)[1]
  high_col <- grep("UCL", names(tab), value = TRUE)[1]
  med <- data.frame(
    group = tab$group,
    median = tab[[med_col]],
    conf.low = tab[[low_col]],
    conf.high = tab[[high_col]],
    stringsAsFactors = FALSE
  )

  raw_counts <- .ekb_survival_raw_counts(object)
  med <- merge(med, raw_counts, by = "group", all.x = TRUE, sort = FALSE)

  landmarks <- sort(unique(as.numeric(landmarks)))
  lm_df <- data.frame()
  if (anyNA(landmarks) || any(!is.finite(landmarks)) || any(landmarks < 0))
    stop("landmarks must be finite, non-negative times.", call. = FALSE)
  if (length(landmarks) && (isTRUE(extend) || any(landmarks <= max(fit$time)))) {
    s <- suppressWarnings(summary(fit, times = landmarks, extend = extend, conf.int = conf.level))
    if (length(s$time)) {
      grp <- if (is.null(s$strata)) rep("Overall", length(s$time)) else sub("^[^=]+=", "", as.character(s$strata))
      lm_df <- data.frame(
        group = grp,
        time = s$time,
        estimate = s$surv,
        conf.low = s$lower,
        conf.high = s$upper,
        n.risk = s$n.risk,
        n.event = s$n.event,
        n.censor = s$n.censor,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    median = .ekb_as_tibble(med),
    landmarks = .ekb_as_tibble(lm_df),
    raw_counts = .ekb_as_tibble(raw_counts),
    weighted = !is.null(object$weights),
    conf.level = conf.level,
    extend = extend
  )
}

.ekb_risk_table_data <- function(object, times) {
  dat <- object$data
  cols <- c(object$time, object$event, object$group)
  cols <- stats::na.omit(cols)
  keep <- rep(TRUE, nrow(dat))
  if (!is.null(object$fit$na.action)) keep[as.integer(object$fit$na.action)] <- FALSE
  dat <- dat[keep, , drop = FALSE]

  if (is.null(object$group)) {
    groups <- "Overall"
    dat$.ekb_group <- "Overall"
  } else {
    dat$.ekb_group <- as.character(dat[[object$group]])
    groups <- levels(droplevels(factor(dat$.ekb_group)))
  }

  out <- do.call(rbind, lapply(groups, function(g) {
    d <- dat[dat$.ekb_group == g, , drop = FALSE]
    data.frame(
      group = g,
      time = times,
      n.risk = vapply(times, function(tt) sum(d[[object$time]] >= tt, na.rm = TRUE), integer(1)),
      stringsAsFactors = FALSE
    )
  }))
  out$group <- factor(out$group, levels = rev(groups))
  out
}

#' Plot survival curves and risk counts
#'
#' Plots an ekb_survival_fit object with optional confidence bands and raw risk counts.
#'
#' @param object A fitted ekbMed object appropriate to this function; see the description.
#' @param xlab Horizontal axis label.
#' @param ylab Vertical axis label.
#' @param xlim Optional numeric pair of time-axis limits.
#' @param break.time.by Positive spacing between time-axis ticks.
#' @param conf.int Whether to display confidence intervals.
#' @param censor Whether to mark censoring times.
#' @param risk.table Whether to add an unweighted number-at-risk panel.
#' @param show.p Whether to display the log-rank p-value.
#' @param legend.title Optional legend title.
#' @param palette Optional vector of group colours.
#' @param linetypes Whether to distinguish groups by line type.
#' @param percent Whether the vertical axis shows percentages.
#' @param text.size Base plot text size in points.
#' @return A ggplot object, or a patchwork composition when risk.table is TRUE.
#' @details Requires ggsurvfit for step ribbons and patchwork when risk.table is TRUE. Weighted p-values require RISCA. The risk-table heading explicitly identifies raw unweighted counts when curves use weights.
#' @export
#' @examples
#' set.seed(21)
#' d <- data.frame(
#'   time = rexp(80, 0.04), event = rbinom(80, 1, 0.7),
#'   treatment = factor(rep(c("A", "B"), 40)), age = rnorm(80, 60, 8),
#'   sex = factor(rep(c("F", "F", "M", "M"), 20))
#' )
#' if (requireNamespace("ggsurvfit", quietly = TRUE)) {
#'   plot_survival(fit_survival(d, "time", "event", "treatment"), risk.table = FALSE)
#' }
plot_survival <- function(object,
                          xlab = "Time (months)",
                          ylab = "Survival probability",
                          xlim = NULL,
                          break.time.by = 12,
                          conf.int = TRUE,
                          censor = TRUE,
                          risk.table = TRUE,
                          show.p = TRUE,
                          legend.title = NULL,
                          palette = NULL,
                          linetypes = TRUE,
                          percent = TRUE,
                          text.size = 11) {
  .ekb_require("ggplot2")
  .ekb_require("broom")
  .ekb_require("scales")
  .ekb_require("ggsurvfit")

  if (!inherits(object, "ekb_survival_fit")) stop("object must be created by fit_survival().", call. = FALSE)

  fit0 <- survival::survfit0(object$fit)
  df <- broom::tidy(fit0)
  if (!"strata" %in% names(df)) df$strata <- "Overall"
  df$strata <- sub("^[^=]+=", "", as.character(df$strata))

  aes_base <- ggplot2::aes(x = .data$time, y = .data$estimate, group = .data$strata,
                           color = .data$strata)
  p <- ggplot2::ggplot(df, aes_base)

  if (isTRUE(conf.int) && all(c("conf.low", "conf.high") %in% names(df))) {
    p <- p + ggsurvfit::stat_stepribbon(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$strata),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    )
  }

  if (isTRUE(linetypes) && length(unique(df$strata)) > 1L) {
    p <- p + ggplot2::geom_step(ggplot2::aes(linetype = .data$strata), linewidth = 0.9)
  } else {
    p <- p + ggplot2::geom_step(linewidth = 0.9)
  }

  if (isTRUE(censor) && "n.censor" %in% names(df)) {
    cens <- df[df$n.censor > 0, , drop = FALSE]
    if (nrow(cens)) {
      p <- p + ggplot2::geom_point(
        data = cens,
        ggplot2::aes(x = .data$time, y = .data$estimate, color = .data$strata),
        inherit.aes = FALSE,
        shape = 3,
        size = 1.8,
        stroke = 0.7,
        show.legend = FALSE
      )
    }
  }

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette) + ggplot2::scale_fill_manual(values = palette)
  }

  y_scale <- if (isTRUE(percent)) {
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = c(0, 0))
  } else {
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  }

  xmax <- if (!is.null(xlim)) max(xlim) else max(df$time, na.rm = TRUE)
  xmin <- if (!is.null(xlim)) min(xlim) else 0
  breaks <- seq(xmin, xmax, by = break.time.by)

  p <- p +
    y_scale +
    ggplot2::scale_x_continuous(breaks = breaks, limits = c(xmin, xmax), expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = ylab, color = legend.title, linetype = legend.title) +
    theme_ekbmed(base_size = text.size) +
    ggplot2::theme(legend.position = "top")

  if (isTRUE(show.p) && !is.null(object$group)) {
    lr <- survival_logrank(
      object$data,
      time = object$time,
      event = object$event,
      group = object$group,
      weights = object$weights
    )
    p_label <- if (lr$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", lr$p.value)
    p <- p + ggplot2::annotate("text", x = xmin + 0.03 * (xmax - xmin), y = 0.08,
                               label = p_label, hjust = 0, size = text.size / 3)
  }

  if (!isTRUE(risk.table)) return(p)
  .ekb_require("patchwork")
  risk <- .ekb_risk_table_data(object, breaks)
  rt <- ggplot2::ggplot(risk, ggplot2::aes(x = .data$time, y = .data$group, label = .data$n.risk)) +
    ggplot2::geom_text(size = text.size / 3) +
    ggplot2::scale_x_continuous(breaks = breaks, limits = c(xmin, xmax), expand = c(0, 0)) +
    ggplot2::labs(
      x = xlab,
      y = NULL,
      title = if (is.null(object$weights)) "Number at risk" else "Number at risk (unweighted)"
    ) +
    theme_ekbmed(base_size = max(9, text.size - 1)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "plain", size = max(9, text.size - 1)),
      legend.position = "none"
    )

  p_no_x <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  p_no_x / rt + patchwork::plot_layout(heights = c(4, 1.25))
}
