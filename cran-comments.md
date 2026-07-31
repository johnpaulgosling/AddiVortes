## Submission of AddiVortes 0.6.9

This resubmission addresses two CRAN check failures from earlier 0.6.8/0.6.9 builds:

1. AddressSanitizer heap-buffer-overflow in `propose_internal()`
   (spherical vignette rebuild on SAN / clang-ASAN builders).
2. Install WARNING from unused variable `-Wunused-variable` in
   `log_acceptance_components()`.

## Test environments

* local Ubuntu 24.04, R 4.6.1
* previous CRAN SAN / clang-ASAN reports for 0.6.8 (ASan abort, now fixed)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

This flag is injected by the system R installation on Ubuntu
(`/usr/lib/R/etc/Makeconf`), not by this package. The package has no
`src/Makevars` and sets no compiler flags of its own. The NOTE does not
appear on CRAN's builders.
