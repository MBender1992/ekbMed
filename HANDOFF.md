# Handoff: ekbMed 0.1.0

Stand: 26. September 2026. Das ZIP enthält den direkt als Repository nutzbaren
Ordner `ekbMed/`, ohne `.git/` und ohne private Patienten- oder Studiendaten.

## Lieferumfang

- Reguläres R-Paket: DESCRIPTION, roxygen2-generierter NAMESPACE, R/ und man/.
- README.md mit simulierten Beispielen für die vereinheitlichte API.
- SPECIFICATIONS.md mit Architektur, statistischen Konventionen und Grenzen.
- NEWS.md mit den minimalen Integrationskorrekturen.
- SOURCE_MANIFEST.csv mit Herkunft und SHA-256 der freigegebenen Ausgangsskripte.
- tests/testthat mit 44 Testblöcken und ausschließlich synthetischen Daten.
- GitHub Actions für Dokumentationsabgleich und R CMD check.
- validation/ mit tatsächlichen Prüfprotokollen und Dependency-Versionen;
  dieser Ordner wird beim R-Paket-Build ausgeschlossen.

## Öffentliche API

46 exportierte Funktionen; zusätzlich registrierte S3-Methode
`print.ekb_clinical_qc`. Interne `.ekb_*`-Helper sind nicht exportiert.

`add_median_survival`, `as_ekb_flextable`, `ate_weights`, `calc_DCR`, `calc_ORR`, `calc_response_rate`, `calc_survival`, `calc_survival_interval`, `check_cox_ph`, `convert_date`, `convert_response`, `cox_global_wald`, `cox_output`, `coxph_meta_analysis`, `diagnose_weighting`, `fit_cox`, `fit_flextable_to_page`, `fit_iptw_survival`, `fit_mi_cox`, `fit_mi_iptw_cox`, `fit_survival`, `fit_weighting`, `format_cox`, `impute_clinical`, `manual_ate_weights`, `plot_balance`, `plot_cox_forest`, `plot_iptw_survival`, `plot_subgroup_forest`, `plot_survival`, `qc_clinical_data`, `subgroup_cox`, `subgroup_mi_cox`, `summarize_survival`, `survival_logrank`, `survival_time`, `tbl_baseline`, `tbl_cox`, `tbl_outcomes`, `tbl_subgroup`, `tbl_survival`, `tbl_treatment_course`, `theme_ekbmed`, `theme_ekbmed_flextable`, `tidy_cox`, `truncate_weights`.

## Abhängigkeiten

Die tatsächlichen Imports und optionalen Suggests stehen in DESCRIPTION;
`validation/dependency-versions.csv` dokumentiert die geprüften Versionen.
Survival, Datenaufbereitung und grundlegende Grafik sind Imports. Gewichtung,
MI, Tabellen/Word-Export, spezielle Grafiken und RISCA sind optionale Backends
in Suggests mit Laufzeitprüfung. Für den vollständigen lokalen Check alle
Suggests installieren; siehe README. Entwicklungswerkzeuge sind keine Imports.

## Ausgeführte Validierung

- R 4.3.3, Linux x86_64 (Ubuntu 24.04).
- `devtools::document("ekbMed")`: erfolgreich; man/ und NAMESPACE erzeugt.
- `devtools::test("ekbMed")`: 44 Testblöcke, 224 bestandene Assertions,
  0 fehlgeschlagene Assertions, 0 Warnungen, 1 übersprungener Test.
- Der Skip betrifft ausschließlich den direkten Vergleich des gewichteten
  Log-rank-Tests mit RISCA: dieses optionale Paket war nicht installiert.
- Cox- und Subgruppen-Forestplots zusätzlich gerendert und visuell geprüft.
- `R CMD build`: erfolgreich.
- `R CMD check --no-manual --no-build-vignettes`: **0 Errors, 0 Warnings,
  1 Note**. Die einzige Note nennt das nicht installierte Suggests-Paket RISCA.
  Hierfür wurde `_R_CHECK_FORCE_SUGGESTS_=false` verwendet. Alle anderen
  Prüfungen einschließlich ausführbarer Beispiele und installierter Tests
  bestanden. Das PDF-Handbuch wurde mangels LaTeX nicht erzeugt.

`devtools::install(..., dependencies = FALSE, upgrade = "never", build = FALSE)`
und ein Smoke-Test über `library(ekbMed)` waren erfolgreich. Dabei meldete die
portable Ausführungsumgebung einmal, dass eine Prozess-Stat-Datei nicht gelesen
werden konnte; Installation und Funktionsaufruf wurden dennoch abgeschlossen.
Diese Umgebungswarnung trat nicht im abgeschlossenen R CMD check auf.

Die Tests vergleichen Cox-Modelle und robuste Standardfehler direkt mit
survival, Gewichte mit WeightIt und MI-Pooling einschließlich robuster
Within-Imputation-Varianz mit mice. Sie prüfen beide MI-Subgruppengewichtungen,
Wiederverwendung der Imputation, Messages, Response-Nenner, Tabellenheader,
Total-Spalten, Missingness, Gewichtungsmetadaten und darstellbare Grafiken.
Die historische PDSeq-Validierung stammt aus dem Projektkontext; sie wurde
hier nicht mit privaten Originaldaten wiederholt.

## Minimale Änderungen gegenüber den freigegebenen Skripten

1. Baseline-Tabellen nehmen die intern ergänzte Gewichtsspalte nicht mehr
   versehentlich in eine implizite Variablenauswahl auf.
2. Baseline-Gruppen als Zeichenwerte erhalten übereinstimmende Reihenfolge
   in Header und statistischem Tabelleninhalt.
3. Flextable-Konvertierung akzeptiert auch gtsummary-Objekte ohne ekb-Metadaten.
4. Outcome-Fußnoten zeigen das angeforderte Konfidenzniveau.
5. Cox-Header verwenden die unterstützte gtsummary-Spalte `conf.low`.
6. Unicode-Escapes, Namespace-Imports, S3-Registrierung und dokumentierte
   Data-mask-Symbole stellen die Paketkompatibilität her.
7. Neue roxygen-Beispiele wurden im Check korrigiert und erneut generiert.

Keine Änderungen an Estimands, Referenzgruppen, statistischen Definitionen,
Missing-Data-Logik, Imputationsarchitektur oder Varianzdefaults.
MIT wurde entsprechend der ausdrücklichen Erlaubnis im Auftrag gewählt;
Autorenname und E-Mail stammen aus der bestehenden Repository-DESCRIPTION.

## Lokaler Workflow

```r
library(devtools)
library(roxygen2)
# Im übergeordneten Ordner des entpackten ekbMed-Verzeichnisses:
package <- "ekbMed"
setwd(package)
document()
setwd("..")
check(package)
install(package)
```

Das Paket enthält keine absoluten lokalen Pfadannahmen. Ein vollständiger
lokaler Standard-Check benötigt alle Suggests sowie die üblichen R-Build-
Werkzeuge (und für das PDF-Handbuch eine LaTeX-Installation).
Das bestehende GitHub-Repository wurde nicht verändert.
