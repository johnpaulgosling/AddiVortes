## Submission of AddiVortes 1.0.2

This is a patch release. `AddiVortes()` gains `eskinWeight`, a positive
weight on the Eskin mismatch used when `cat.onehot = FALSE`. The default
weight is 1, so existing Eskin fits are unchanged. Setting the weight to
`d^2 / 2` makes every categorical mismatch contribute 1. The training
level counts and weights are stored on the fit and reused by `predict()`.

## Test environments

* GitHub Actions: macOS-latest (R-release), Windows-latest (R-release),
  Ubuntu-latest (R-devel, R-release, R-oldrel)
* local Ubuntu, R 4.3.3

## R CMD check results

0 errors | 0 warnings | 0 notes on CRAN-like builders.

On some Ubuntu and macOS installations, `R CMD check --as-cran` reports:

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Werror=format-security’ ‘-Wformat’ ‘-Wp,-D_FORTIFY_SOURCE=3’
    ‘-Wp,-D_GLIBCXX_ASSERTIONS’ ‘-march=x86-64-v3’ ‘-mpclmul’

Those flags come from the system R `Makeconf` (for example
`/usr/lib/R/etc/Makeconf` on Ubuntu, or the macOS R binary used on GitHub
Actions). This package has no `src/Makevars` and sets no compiler flags of
its own. The NOTE does not appear on CRAN's builders. GitHub Actions sets
`_R_CHECK_COMPILATION_FLAGS_KNOWN_` so the same NOTE is not treated as a
check issue on this repository's runners.
