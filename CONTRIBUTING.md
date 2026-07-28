# Contributing to multihostSIR

Thanks for your interest in the project.

## Reporting bugs

Open an issue at https://github.com/kramera3/multihostSIR/issues with a minimal
reproducible example and the output of `sessionInfo()`. A `reprex::reprex()` is
ideal.

## Suggesting features

Open an issue describing the ecological or epidemiological question the feature
would answer. That context helps more than an API sketch.

## Pull requests

1. Fork and create a branch from `main`.
2. Make your change. Add or update tests in `tests/testthat/`.
3. Add roxygen documentation (`#'` comments) for any new or changed function.
   If the change is user-facing, consider whether a vignette in `vignettes/`
   should mention it too.
4. Run:

   ```r
   devtools::document()
   devtools::test()
   devtools::check()   # builds the vignettes; needs pandoc
   ```

   All three must be clean before you open the PR.
5. Add a bullet to `NEWS.md` under a "development version" heading.
6. Open the pull request describing what changed and why.

## Style

* Follow the tidyverse style guide: https://style.tidyverse.org/
* Keep lines under 80 characters.
* Use `snake_case` for functions and variables.
* Comments explain *why*, not *what*. Assume the reader can read R.
* New exported functions need examples that run in under 5 seconds.

## Scope

The package models multi-host SIR dynamics with allometric parameterisation.
Additions that broaden it substantially (different compartmental structures,
spatial dynamics, stochastic individual-based models) are better as separate
packages that build on this one. Open an issue to discuss first.
