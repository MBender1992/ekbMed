# ekbMed technical specification

Version 0.1.0.

## Philosophy and scope

ekbMed provides reusable tools for clinical survival analysis and reporting.
Its API separates data preparation, propensity-score weighting, multiple
imputation, outcome modelling and presentation. Analysis choices are explicit:
researchers specify covariates, reference groups, estimands and missing-data
strategies according to their study design and statistical analysis plan.

Cohort eligibility, endpoint construction, censoring rules and data linkage
belong in the analysis project. The package does not choose these definitions
or perform automatic model selection. No clinical datasets are distributed.

This document describes object contracts, statistical conventions and developer
expectations. See [README.md](README.md) for worked examples, function help for
arguments and return values, and [NEWS.md](NEWS.md) for version changes.

## Public API and object contracts

| Producer | Class / structure | Consumers |
|---|---|---|
| fit_survival | ekb_survival_fit: fit, data, endpoint names, numeric weights | summarize_survival, plot_survival |
| fit_cox | ekb_cox_fit: coxph fit, input, covariates, weights/weighting | tidy_cox, check_cox_ph, tables, forest plots |
| fit_weighting | ekb_weighting: WeightIt object, data, weights, propensity_score, estimand | diagnose_weighting and weights arguments |
| subgroup_cox | ekb_subgroup_cox: results, optional interaction_tests, reference/comparator | tbl_subgroup, plot_subgroup_forest |
| impute_clinical | ekb_imputation: mids, input data, predictor matrix, methods and seed | fit_mi_cox, subgroup_mi_cox |
| fit_mi_cox | ekb_mi_cox: original imputation, fits, pooled, pooled_summary, diagnostics | tbl_cox, plot_cox_forest |
| subgroup_mi_cox | ekb_subgroup_mi_cox: results, original imputation, global weighting, diagnostic_summary | tbl_subgroup, plot_subgroup_forest |
| fit_iptw_survival | ekb_iptw_survival: adjustedCurves fit and metadata | plot_iptw_survival |
| qc_clinical_data | ekb_clinical_qc: overview, variable and issue summaries | print method, analyst review |
| tbl_baseline / tbl_cox / tbl_subgroup / tbl_outcomes / tbl_treatment_course | gtsummary plus ekb attributes | as_ekb_flextable |
| tbl_survival | tibble | standard tibble/flextable consumers |

Functions are documented individually under `man/`; the generated NAMESPACE is
the definitive export list. Internal `.ekb_*` functions and `%||%` are unexported.
The QC print method is registered as an S3 method.

## Weight object handling

The shared resolver accepts NULL, a numeric vector, a column name, or an
`ekb_weighting` object. It validates length and positive finite nonmissing
weights. NA weights remain subject to the consumer's missingness policy;
baseline tables require complete weights. Matching is positional, not by ID.
Do not silently join/reorder observations between weighting and outcome models.

Ordinary-dataset weighting is upstream and reusable. Weighted Cox can contain
additional outcome covariates, which changes the interpretation relative to a
weighted-only marginal model. No doubly-robust claim is made. Truncation is an
explicit helper and never automatic.

## Imputation architecture

`impute_clinical` constructs mice once with explicit targets. Targets can predict
each other; additional predictors and auxiliary variables are explicit.
Treatment and outcome fields cannot be targets. The survival option includes
event and Nelson–Aalen hazard, not raw survival time by default. Unimputed
predictors with missing values are rejected.

Downstream functions accept ekb_imputation or mice::mids and call complete for
each dataset. No downstream function calls mice again. With PS covariates,
WeightIt is fitted within every completed dataset and each fit receives its own
weights. Model estimates are pooled with mice::pool. Diagnostic summaries expose
imputation-specific balance, ESS and positivity flags.

A completed imputation is not a complete-case dataset in the analytic sense.
Correctness of the imputation model and MAR assumptions remain the analyst's
responsibility. Seeds and model metadata are retained for reproducibility.

## Subgroup weighting

Global scope estimates full-cohort weights once per imputation and subsets them
for each level. Within-subgroup scope estimates weights separately within each
level and imputation after removing the subgroup and constant predictors. A
varying predictor must remain. The overall row always uses full-cohort weights.
This implementation calculates global diagnostic fits even in local scope.

Outcome adjustment removes the subgroup variable and locally constant
covariates. Single-dataset sparse levels return a note and unavailable estimates;
MI levels must be represented with two treatment groups in every imputation.
MI subgroup counts are the mean over imputations; displayed integer N is rounded.
The API does not pool interaction tests. Default subgroup p-values are tests of
treatment effects within levels.

## Missingness, endpoints and variance

Event coding must be 0/1 or logical. Times must be finite and non-negative;
counting-process start must precede stop. Cox missingness follows its backend
and session na.action; survival defaults explicitly to na.omit. PS fitting
requires complete baseline predictors and treatment. No missing-data rule is
silently replaced by imputation or a newly chosen estimator.

`robust = NULL` delegates to coxph; non-integer weights generally trigger sandwich
variance. `robust = TRUE` requests it explicitly, including integer weights.
`tidy_cox` uses robust.se when available. Mice >= 3.13.2 pools robust.se from broom
when present. Tests compare the pooled within-imputation variance against direct
sandwich variances. Weights are treated as supplied case weights; no extra
propensity-estimation variance adjustment is introduced.

Date parsing retains dmy-first priority and explicit partial-date rules.
Inclusive intervals default to +1 day; no RFS/TTNT/TOT rule is defined here.
Response conversion provides configurable NC/NED/MR conventions.
ORR and DCR standalone helpers exclude non-evaluable response codes by default.

## Reporting conventions

Survival plots use raw observation risk counts, clearly labelled as unweighted
when curves are weighted. Landmark backend risk counts can be weighted.
Medians and landmarks are not extrapolated beyond follow-up by default.
Weighted log-rank uses RISCA and requires two groups; unweighted uses survdiff.

Baseline tables include Total by default; weighted headers explicitly use
Weighted N and define it as sum of weights. Outcome tables keep observed BOR and
binary response descriptions unweighted while optionally weighting survival.
Their observed-N headers and weighting footnotes express this distinction.
Response binary columns are created explicitly by the analyst, preserving NA.
Month-labelled tables assume follow-up in months.

Treatment-course tables treat blanks as missing, support explicit level order,
and summarize toxicity details within toxicity-positive patients by default.
Flextables are editable, use black-and-white styling, bold variable blocks,
indented levels and more deeply indented landmarks. `tbl_survival` returns a
tibble; it is not routed through as_ekb_flextable directly.
Forest plots are ggplot objects with level-effect p-values by default.

## Messages and dependencies

Single-dataset Cox/subgroup calls announce ekb_weighting input once. MI calls
announce within-imputation IPTW once, including scope for subgroup models.
Internal Cox loops receive numeric weights. WeightIt informational messages in
MI loops are suppressed, but its warnings/errors are retained. Outcomes and
weighted baseline tables announce their interpretation once.

Core namespaces are Imports; heavier feature backends are Suggests and checked
at the point of use. All package functions avoid library/require attachment.
NAMESPACE and Rd documentation are generated by roxygen2, never maintained by
hand. Tests contain no network operations and skip clearly when optional
features are unavailable. CI installs dependencies before checking.

## Compatibility

Historical emR wrappers delegate without their own statistical implementation:
ate_weights, survival_time, add_median_survival, cox_output,
coxph_meta_analysis, fit_mi_iptw_cox; calc_survival retains its argument spelling.
The combined legacy MI-IPTW wrapper warns and calls the two current entry points.
Automatic backward selection is not implemented: cox_output warns for requested
non-full models and fits the full model. Compatibility wrappers are documented separately from the recommended API.

## Deliberate exclusions

The package does not automatically select models, trim weights, change
reference groups, impute outcomes or orchestrate study-specific reports.
Subgroup effect estimates do not establish treatment-effect heterogeneity;
interaction testing requires a separate analysis decision. Weighting and
imputation do not establish causal identification, the validity of missing-at-
random assumptions or double robustness by themselves.

## Testing and interpretation

The test suite uses synthetic data and checks results against direct calls to
survival, WeightIt and mice. Coverage includes unweighted and weighted Cox
models, robust variance pooling, reusable imputations, subgroup weighting
scopes, response summaries, tables, plots, input validation and messages.
Tests require no network access and report skips when an optional backend is
unavailable. Full release checks should install all suggested dependencies.

These checks assess software behaviour and agreement with the specified
backends. They do not establish the suitability of an analysis for a particular
study. Researchers must assess model assumptions, covariate balance, overlap,
imputation diagnostics and the clinical interpretation of their endpoints.

## Maintenance conventions

- Generate NAMESPACE and man/ from roxygen2 comments in R/.
- Accompany behavioural changes with focused tests and a NEWS entry.
- Document changes to defaults, estimands, missingness or variance handling
  explicitly; do not introduce them as incidental implementation changes.
- Keep examples reproducible and independent of private data or local paths.
- Record package versions and analysis settings alongside study results.
