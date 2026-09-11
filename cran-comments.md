## Submission of AddiVortes 1.0.1

This is a patch release. It hides the internal `cellIndices()` helper,
updates the default MCMC settings to 2000 iterations with 500 burn-in,
adds Leo Thomson as an author, and addresses the clang
`-Wabsolute-value` WARNING reported on CRAN's
`r-devel-linux-x86_64-debian-clang` builder.

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
