# Copilot Instructions for AddiVortes

## Project Overview

AddiVortes is an R package implementing **Bayesian Additive Voronoi Tessellations** (BAVT), a non-parametric regression model using Voronoi tessellations as an alternative to BART (Bayesian Additive Regression Trees). The package supports Euclidean and spherical metrics, categorical covariates, and posterior inference via MCMC.

Reference: Stone, A. and Gosling, J.P. (2025). *AddiVortes: (Bayesian) Additive Voronoi Tessellations*. Journal of Computational and Graphical Statistics. <doi:10.1080/10618600.2024.2414104>

## Repository Structure

- `R/` — R source files. Files prefixed with `0_` are internal helpers; `AddiVortes.R` is the main fitting function; `AddiVortesFit.R` defines the S3 class and its methods.
- `src/` — C++ source (`addi_vortes_code.cpp`) for performance-critical routines, registered via `src/init.c`.
- `tests/testthat/` — Unit tests using the `testthat` (edition 3) framework.
- `vignettes/` — R Markdown vignettes (`introduction.Rmd`, `prediction.Rmd`).
- `man/` — Auto-generated Roxygen2 documentation (do not edit by hand).
- `.github/workflows/` — CI workflow (`check-install.yml`) that installs and checks the package on Ubuntu with both GCC and Clang.

## Coding Conventions

### R Code
- Follow the tidyverse style guide for R code (spaces around operators, 2-space indentation).
- Use `camelCase` for variable and function names (e.g., `xScaled`, `proposeTessellation`).
- All exported functions must have complete Roxygen2 documentation (`@title`, `@description`, `@param`, `@return`, `@examples`).
- Internal (non-exported) helper functions have names ending in `_internal` (e.g., `scaleData_internal`, `encodeCategories_internal`).
- Use `\donttest{}` in `@examples` for slow-running examples.
- Do not use `library()` or `require()` inside package functions; use `::` for external package calls or declare in `NAMESPACE`.

### C++ Code
- Performance-critical Voronoi/kNN computations are implemented in `src/addi_vortes_code.cpp`.
- Follow existing patterns when adding or modifying C++ routines; register new `.Call` entry points in `src/init.c`.

### Tests
- Tests use `testthat` edition 3; place new test files in `tests/testthat/` with names matching `test-*.R` (hyphenated, e.g., `test-my-feature.R`).
- Always call `set.seed()` before any stochastic operations in tests to ensure reproducibility.
- Use `showProgress = FALSE` when calling `AddiVortes()` inside tests to suppress output.
- Prefer `expect_equal(round(..., 3), expected_value)` for numeric comparisons involving MCMC output.

## Key Algorithmic Details

- **MCMC loop**: The main backfitting algorithm is in `R/AddiVortes.R`. Each iteration proposes a new tessellation via `R/0_ProposeTessellation.R` and updates cell means via `R/0_SampleMuValues.R`.
- **Metrics**: `metric = 0` → Euclidean distance; `metric = 1` → Spherical distance. Per-dimension metric vectors are supported.
- **Categorical covariates**: Encoded via `encodeCategories_internal()` in `R/0_CategoricalEncoding.R`. Categories become one-hot binary columns scaled by `catScaling`. The encoding metadata is stored in the `catEncoding` field of the fitted object.
- **Scaling**: Continuous covariates are shift-scaled to `[0, 1]` via `scaleData_internal()` in `R/0_ShiftScale.R`. Spherical dimensions are handled with periodic wrapping. Binary (categorical) columns bypass standard scaling.
- **Parallelism**: Parallel chains are supported via the `parallel` package.

## Building and Testing

- **Install dependencies**: `Rscript -e 'remotes::install_deps(dependencies = TRUE)'`
- **Run tests**: `Rscript -e 'testthat::test_package("AddiVortes")'`
- **Build documentation**: `Rscript -e 'roxygen2::roxygenise()'`
- **Check package**: `R CMD check .`
- R (>= 3.5.0) is required; `parallel` (>= 4.0.0) and `pbapply` (>= 1.6) are runtime dependencies.

## Pull Request Guidelines

- Keep changes focused and minimal; avoid unrelated refactoring.
- Ensure `R CMD check` passes with no errors or warnings before opening a PR.
- Add or update tests for any new functionality.
- Update Roxygen2 documentation and re-run `roxygenise()` when changing function signatures.
- Target the `main` branch for all PRs.
