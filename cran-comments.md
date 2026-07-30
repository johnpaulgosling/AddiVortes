## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’

This flag is injected by the system R installation on Ubuntu
(`/usr/lib/R/etc/Makeconf`), not by this package. The package has no
`src/Makevars` and sets no compiler flags of its own. The NOTE does not
appear on CRAN's builders.
