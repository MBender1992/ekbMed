#' Reproducible clinical survival analysis and reporting
#'
#' A unified API for survival models, propensity-score weighting, multiple
#' imputation and clinical reporting. Start with [fit_cox()], [fit_survival()],
#' [fit_weighting()] or [impute_clinical()]. Use [fit_mi_cox()] and
#' [subgroup_mi_cox()] to reuse an upstream imputation.
#'
#' @details Covariates and endpoints must be prespecified by the analyst.
#'   The package does not select models or define study-specific eligibility.
#'   Inspect propensity-score overlap and balance before interpreting weighted
#'   results. Conditional, additionally adjusted hazard ratios need not equal
#'   marginal weighted-only hazard ratios.
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom graphics plot
#' @importFrom utils head tail
#' @keywords package
"_PACKAGE"

# Names evaluated in gtsummary data masks and tidyselect expressions.
utils::globalVariables(c("p.value", "label", "is_heading", "n", "p_value",
                         "subgroup", "hr_ci"))
