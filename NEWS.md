# ekbMed 0.1.0

## Package integration

- Package structure, generated roxygen2 help/NAMESPACE, README, specifications,
  source manifest, synthetic-data testthat suite and R CMD check CI.
- Latest approved API refactor and subsequent IPTW message update retained.
- Publication table module uses the approved v6 source, including Total columns,
  weighted headers, outcome survival weighting and ordered toxicity categories.
- Historical wrappers retained as delegating compatibility functions.
- No changes to estimands, references, imputation models, response definitions,
  survival endpoint definitions or variance defaults.
- MIT license adopted under the owner's package-conversion instruction; author
  identity retained from the repository DESCRIPTION.

## Minimal integration fixes

- Weighted baseline tables with include = NULL select original data columns
  before adding the internal weight column. The internal weight is no longer
  eligible for accidental display as a baseline variable.
- The flextable renderer handles gtsummary objects without ekb_table_type
  metadata without a length-zero logical condition.
- Non-ASCII characters in executable R strings use portable Unicode escapes;
  displayed labels retain the original typography.
- Register the clinical-QC print method and import symbols used within the
  package namespace; declare documented gtsummary data-mask names for checks.

See HANDOFF.md for the actual validation outcome and any further integration
fixes identified during execution.

- Character grouping columns in baseline tables now preserve the same order
  in the statistical body and column headers, preventing mislabelled groups.
- Outcome survival footnotes display the requested confidence level instead of
  a hard-coded 95 percent. Estimates themselves are unchanged.
- Cox table headers use the supported conf.low column rather than gtsummary's
  deprecated ci column.
