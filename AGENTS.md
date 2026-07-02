# AGENTS.md

## Cursor Cloud specific instructions

`AddiVortes` is a standalone R package (R + C++20 via `src/addi_vortes_code.cpp`); there are no servers, databases, ports, or secrets. "Running the application" means loading the library and fitting/predicting a model. The startup update script already refreshes R dependencies (from `DESCRIPTION`) and reinstalls the package with `R CMD INSTALL .`, so the compiled `AddiVortes.so` is available on session start.

Non-obvious caveats:

- System tools (R >= 4.0.0, a C++20 compiler `g++`, and `pandoc`) are provided by the VM snapshot, not the update script. If R is missing, install via `sudo apt-get install -y r-base-dev` (and `pandoc` for vignette/README checks).
- R packages install fast from the Posit binary mirror for Ubuntu noble: set `options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))` plus a matching `HTTPUserAgent`, otherwise CRAN serves slow source builds.
- Run the test suite with `Rscript -e 'testthat::test_local(".")'` (187 tests). Do NOT use `testthat::test_check("AddiVortes")` against the installed package: the package is installed without tests (no `--install-tests`), so it reports "No test files found".
- Lint config lives in `.lintr.R` (not the usual `.lintr`) and defines a `linters` object. Run it with `Rscript -e 'library(lintr); source(".lintr.R"); print(lint_package(linters = linters))'`. There are pre-existing style lints, mostly in `vignettes/`.
- Full CI-equivalent check: `R CMD build . --no-build-vignettes` then `R CMD check --no-manual --as-cran --ignore-vignettes <tarball>`. Expect only benign NOTEs (non-portable default R compile flag `-mno-omit-leaf-frame-pointer`, network-time, CRAN-incoming). CI errors only on warnings, not notes.
- `data/Boston.rda` uses generic column names: predictors `x1`..`x13` and response `y` (there is no `medv` column).
- Core API: `AddiVortes(y, x, ...)` returns an object of class `AddiVortes`; predict with `predict(fit, newx, showProgress = FALSE)`. Set `showProgress = FALSE` for non-interactive runs.
