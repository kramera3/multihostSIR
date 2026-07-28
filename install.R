# ============================================================================
# multihostSIR bootstrap installer
# ============================================================================
#
# Checks your R version, installs any missing dependencies, then installs
# multihostSIR from GitHub and verifies that it works.
#
# Run it straight from the web:
#
#   source("https://raw.githubusercontent.com/kramera3/multihostSIR/main/install.R")
#
# Or clone the repository and run:
#
#   source("install.R")
#
# This script deliberately uses only base R. It installs the package without
# building the vignettes, so it needs neither devtools nor pandoc and stays
# fast. The rendered vignettes are on GitHub and CRAN; to build them locally,
# install with remotes::install_github(..., build_vignettes = TRUE) instead.
# ============================================================================

local({

  REPO      <- "kramera3/multihostSIR"
  MIN_R     <- "4.1.0"
  CRAN      <- "https://cloud.r-project.org"

  # package = minimum version ("" means any version)
  REQUIRED  <- c(deSolve = "1.30", ggplot2 = "3.4.0",
                 scales  = "1.2.0", rlang   = "1.0.0")
  OPTIONAL  <- c(openxlsx = "", testthat = "3.0.0")

  say <- function(...) cat(sprintf(...), "\n", sep = "")
  rule <- function() cat(strrep("-", 68), "\n")

  rule()
  say("multihostSIR installer")
  rule()

  # -- 1. R version --------------------------------------------------------
  current_r <- paste0(R.version$major, ".", R.version$minor)
  if (utils::compareVersion(current_r, MIN_R) < 0) {
    stop(sprintf(paste0("R %s or newer is required; you are running %s.\n",
                        "Please upgrade R from https://cloud.r-project.org"),
                 MIN_R, current_r), call. = FALSE)
  }
  say("R version    : %s (ok, needs >= %s)", current_r, MIN_R)
  say("Library path : %s", .libPaths()[1])

  if (!dir.exists(.libPaths()[1]) ||
      file.access(.libPaths()[1], mode = 2) != 0) {
    say("")
    say("WARNING: your primary library is not writable.")
    say("R will normally offer to create a personal library instead.")
  }

  # -- 2. Work out what is missing -----------------------------------------
  needs_install <- function(pkgs) {
    out <- character(0)
    for (p in names(pkgs)) {
      have <- requireNamespace(p, quietly = TRUE)
      if (!have) {
        out <- c(out, p)
        next
      }
      if (nzchar(pkgs[[p]])) {
        v <- as.character(utils::packageVersion(p))
        if (utils::compareVersion(v, pkgs[[p]]) < 0) out <- c(out, p)
      }
    }
    out
  }

  report <- function(pkgs, label) {
    say("")
    say("%s dependencies:", label)
    for (p in names(pkgs)) {
      have <- requireNamespace(p, quietly = TRUE)
      v <- if (have) as.character(utils::packageVersion(p)) else NA
      need <- if (nzchar(pkgs[[p]])) paste0(">= ", pkgs[[p]]) else "any"
      status <- if (!have) "MISSING" else
        if (nzchar(pkgs[[p]]) &&
            utils::compareVersion(v, pkgs[[p]]) < 0) "OUTDATED" else "ok"
      say("  %-10s %-10s %-12s %s", p, need,
          if (is.na(v)) "-" else v, status)
    }
  }

  report(REQUIRED, "Required")
  report(OPTIONAL, "Optional")

  # -- 3. Install what is missing ------------------------------------------
  missing_required <- needs_install(REQUIRED)
  missing_optional <- needs_install(OPTIONAL)
  missing_tools    <- needs_install(c(remotes = ""))

  to_install <- unique(c(missing_required, missing_optional, missing_tools))

  if (length(to_install) > 0) {
    say("")
    say("Installing: %s", paste(to_install, collapse = ", "))
    rule()
    utils::install.packages(to_install, repos = CRAN)
    rule()

    still_missing <- needs_install(REQUIRED)
    if (length(still_missing) > 0) {
      stop(sprintf(paste0("Could not install required package(s): %s\n",
                          "Check the log above. On Linux you may need system ",
                          "libraries; on macOS you may need Xcode command ",
                          "line tools (xcode-select --install)."),
                   paste(still_missing, collapse = ", ")), call. = FALSE)
    }
  } else {
    say("")
    say("All dependencies already satisfied.")
  }

  # -- 4. Install multihostSIR ---------------------------------------------
  say("")
  say("Installing multihostSIR from github.com/%s ...", REPO)
  rule()

  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("The 'remotes' package is required to install from GitHub.",
         call. = FALSE)
  }

  # build_vignettes = FALSE keeps the install fast and avoids needing pandoc.
  # The rendered vignettes are available on GitHub and CRAN; pass TRUE here if
  # you want them built locally and you have pandoc installed.
  remotes::install_github(REPO, upgrade = "never", build_vignettes = FALSE)

  rule()

  # -- 5. Verify -----------------------------------------------------------
  say("")
  say("Verifying the installation ...")

  ok <- tryCatch({
    loadNamespace("multihostSIR")
    res <- multihostSIR::run_multihost_sir(
      richness_levels = 3:4, n_iterations = 3,
      duration = 60, seed = 1, warmup = FALSE)
    isTRUE(nrow(res$results) == 6L) && all(is.finite(res$results$R0))
  }, error = function(e) {
    say("FAILED: %s", conditionMessage(e))
    FALSE
  })

  rule()
  if (isTRUE(ok)) {
    say("multihostSIR %s installed and working.",
        as.character(utils::packageVersion("multihostSIR")))
    say("")
    say("Next steps:")
    say("  library(multihostSIR)")
    say("  check_dependencies()")
    say("  demo(getting_started, package = \"multihostSIR\")")
    say("  ?run_multihost_sir")
  } else {
    say("Installation completed but the smoke test did not pass.")
    say("Please open an issue at")
    say("  https://github.com/%s/issues", REPO)
    say("and include the output of sessionInfo().")
  }
  rule()

  invisible(ok)
})
