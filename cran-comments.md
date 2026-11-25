## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release, so the "New submission" note is expected.

## Test environments

* local: Windows 11, R 4.6.0
* GitHub Actions: windows-latest, macos-latest, ubuntu-latest (R release)
* win-builder: R devel and R release

## Notes

* Examples, tests and vignettes use small iteration counts and at most two
  cores, so they run quickly and within CRAN limits.
* Excel support is optional: writing Excel is guarded behind
  `requireNamespace("openxlsx")` and reading Excel behind
  `requireNamespace("readxl")`. Both packages are in Suggests.

## Downstream dependencies

* There are currently no downstream dependencies.
