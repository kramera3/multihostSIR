#' Dependencies used by multihostSIR
#'
#' @return A list with `required` and `optional` named character vectors, where
#'   names are package names and values are the minimum version (or `""`).
#' @noRd
dependency_table <- function() {
  list(
    required = c(deSolve = "1.30",
                 ggplot2 = "3.4.0",
                 scales  = "1.2.0",
                 rlang   = "1.0.0"),
    optional = c(openxlsx = "",
                 testthat = "3.0.0",
                 remotes  = "")
  )
}

#' What each optional package is for
#' @noRd
optional_purpose <- function() {
  c(openxlsx = "Excel workbook output in export_results()",
    testthat = "running the package test suite",
    remotes  = "installing the package from GitHub")
}

#' Check that dependencies are installed and current
#'
#' Reports which of the packages needed by **multihostSIR** are installed, which
#' are missing, and which are older than the minimum version the package
#' expects. Required packages are normally guaranteed by R at install time;
#' this function is most useful for diagnosing a partially broken library, or
#' for checking the optional packages that unlock extra features.
#'
#' @param include_optional Include optional packages in the report.
#' @param quiet Suppress the printed report and return the table silently.
#'
#' @return Invisibly, a data frame with one row per package and columns
#'   `package`, `type`, `needed`, `installed`, `version` and `ok`.
#'
#' @examples
#' check_dependencies()
#'
#' # Just the machine-readable table
#' deps <- check_dependencies(quiet = TRUE)
#' deps[!deps$ok, ]
#'
#' @seealso [install_dependencies()]
#' @export
check_dependencies <- function(include_optional = TRUE, quiet = FALSE) {
  check_flag(include_optional, "include_optional")
  check_flag(quiet, "quiet")

  deps <- dependency_table()
  pkgs <- deps$required
  type <- rep("required", length(pkgs))
  if (include_optional) {
    pkgs <- c(pkgs, deps$optional)
    type <- c(type, rep("optional", length(deps$optional)))
  }

  installed <- vapply(names(pkgs),
                      function(p) requireNamespace(p, quietly = TRUE),
                      logical(1L))

  version <- vapply(seq_along(pkgs), function(i) {
    if (!installed[i]) return(NA_character_)
    as.character(utils::packageVersion(names(pkgs)[i]))
  }, character(1L))

  ok <- vapply(seq_along(pkgs), function(i) {
    if (!installed[i]) return(FALSE)
    if (!nzchar(pkgs[i])) return(TRUE)
    utils::compareVersion(version[i], pkgs[i]) >= 0
  }, logical(1L))

  out <- data.frame(
    package   = names(pkgs),
    type      = type,
    needed    = ifelse(nzchar(pkgs), paste0(">= ", pkgs), "any"),
    installed = installed,
    version   = version,
    ok        = ok,
    stringsAsFactors = FALSE,
    row.names = NULL)

  if (!quiet) {
    print_dependency_report(out)
  }
  invisible(out)
}

#' @noRd
print_dependency_report <- function(out) {
  mark <- ifelse(out$ok, "ok  ", ifelse(out$installed, "old ", "MISS"))
  cat("multihostSIR dependency check\n\n")
  for (i in seq_len(nrow(out))) {
    cat(sprintf("  [%s] %-10s %-10s %s\n", mark[i], out$package[i],
                out$needed[i],
                if (is.na(out$version[i])) "not installed" else
                  paste("installed", out$version[i])))
  }

  bad_req <- out[!out$ok & out$type == "required", ]
  bad_opt <- out[!out$ok & out$type == "optional", ]

  cat("\n")
  if (nrow(bad_req) == 0L && nrow(bad_opt) == 0L) {
    cat("All dependencies satisfied.\n")
    return(invisible(NULL))
  }
  if (nrow(bad_req) > 0L) {
    cat(sprintf("Missing or outdated REQUIRED: %s\n",
                paste(bad_req$package, collapse = ", ")))
  }
  if (nrow(bad_opt) > 0L) {
    purpose <- optional_purpose()
    cat("Missing or outdated OPTIONAL:\n")
    for (p in bad_opt$package) {
      cat(sprintf("  %-10s needed for %s\n", p,
                  if (p %in% names(purpose)) purpose[[p]] else "extra features"))
    }
  }
  cat("\nRun install_dependencies() to install them.\n")
  invisible(NULL)
}

#' Install missing dependencies
#'
#' Installs any package that [check_dependencies()] reports as missing or older
#' than the minimum version required. Packages that are already current are left
#' untouched.
#'
#' @details
#' Nothing is installed without a decision. In an interactive session you are
#' asked to confirm; in a non-interactive session you must pass `ask = FALSE`
#' explicitly, which makes the function safe to call from scripts that are meant
#' to install things and inert everywhere else.
#'
#' @param include_optional Also install the optional packages.
#' @param ask Prompt before installing. Defaults to [interactive()].
#' @param lib Library directory to install into. Defaults to the first element
#'   of [.libPaths()].
#' @param repos CRAN mirror to install from.
#' @param ... Passed to [utils::install.packages()].
#'
#' @return Invisibly, a character vector of the packages that were installed.
#'
#' @examples
#' \dontrun{
#' # Interactive: prompts before doing anything
#' install_dependencies()
#'
#' # Scripted: required packages only, no prompt
#' install_dependencies(include_optional = FALSE, ask = FALSE)
#' }
#'
#' @seealso [check_dependencies()]
#' @export
install_dependencies <- function(include_optional = TRUE,
                                 ask = interactive(),
                                 lib = .libPaths()[1L],
                                 repos = getOption("repos"),
                                 ...) {
  check_flag(include_optional, "include_optional")
  check_flag(ask, "ask")

  status <- check_dependencies(include_optional = include_optional,
                               quiet = TRUE)
  todo <- status$package[!status$ok]

  if (length(todo) == 0L) {
    message("All dependencies are already satisfied; nothing to install.")
    return(invisible(character(0)))
  }

  if (is.null(repos) || !nzchar(repos[[1L]]) ||
      identical(unname(repos[[1L]]), "@CRAN@")) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }

  message("The following packages will be installed into ", lib, ":\n  ",
          paste(todo, collapse = ", "))

  if (ask) {
    if (!interactive()) {
      stop(paste0("`ask = TRUE` in a non-interactive session. Pass ",
                  "`ask = FALSE` to install without prompting."),
           call. = FALSE)
    }
    answer <- readline("Proceed? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("Cancelled; nothing was installed.")
      return(invisible(character(0)))
    }
  }

  utils::install.packages(todo, lib = lib, repos = repos, ...)

  after <- check_dependencies(include_optional = include_optional,
                              quiet = TRUE)
  still_bad <- after$package[!after$ok]
  if (length(still_bad) > 0L) {
    warning("These packages could not be installed: ",
            paste(still_bad, collapse = ", "),
            ". Check the install log above for the reason.", call. = FALSE)
  } else {
    message("All dependencies satisfied.")
  }

  invisible(setdiff(todo, still_bad))
}
