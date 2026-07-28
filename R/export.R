#' Export results to disk
#'
#' Writes the tables from a run to CSV, and optionally to a single Excel
#' workbook. Accepts either a replicate study from [run_multihost_sir()] or a
#' single deterministic run from [run_fixed_community()]; the set of tables
#' written matches the object.
#'
#' @details
#' Excel sheets are capped at 1,048,576 rows. When a table exceeds that
#' limit, the workbook records a pointer note instead of the full table; the
#' complete data are always written to CSV regardless.
#'
#' Writing Excel requires the \pkg{openxlsx} package, which is a suggested
#' rather than a required dependency.
#'
#' @param x A `"multihost_sir"` object from [run_multihost_sir()], or a
#'   `"multihost_fixed"` object from [run_fixed_community()].
#' @param dir Directory to write into. Created if it does not exist.
#' @param formats Which formats to write. Any of `"csv"` and `"xlsx"`.
#' @param prefix File-name prefix for the written files.
#' @param overwrite Overwrite an existing workbook of the same name.
#'
#' @return Invisibly, a character vector of the paths written.
#'
#' @examples
#' res <- run_multihost_sir(3:4, n_iterations = 3, duration = 100, seed = 1)
#' out <- file.path(tempdir(), "multihost_demo")
#' export_results(res, dir = out, formats = "csv")
#' list.files(out)
#'
#' @export
export_results <- function(x, dir = "multihostSIR_outputs",
                           formats = c("csv", "xlsx"),
                           prefix = "multihostSIR",
                           overwrite = TRUE) {
  formats <- match.arg(formats, choices = c("csv", "xlsx"), several.ok = TRUE)
  if (!is.character(dir) || length(dir) != 1L) {
    stop("`dir` must be a single directory path.", call. = FALSE)
  }

  if (inherits(x, "multihost_sir")) {
    tables <- list(
      results       = x$results,
      trajectories  = x$trajectories,
      timing        = x$timing,
      summary       = summary(x)$by_richness,
      predator      = predator_summary(x),
      competence    = competence_summary(x))
  } else if (inherits(x, "multihost_fixed")) {
    tables <- list(
      results        = x$results,
      species        = x$species,
      trajectory     = x$trajectory,
      species_curves = x$species_curves)
  } else {
    stop(paste0("`x` must be created by run_multihost_sir() or ",
                "run_fixed_community()."), call. = FALSE)
  }

  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  written <- character(0)

  if ("csv" %in% formats) {
    for (nm in names(tables)) {
      path <- file.path(dir, sprintf("%s_%s.csv", prefix, nm))
      write.csv(tables[[nm]], path, row.names = FALSE)
      written <- c(written, path)
    }
  }

  if ("xlsx" %in% formats) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("Package 'openxlsx' is not installed; skipping the Excel export.",
              call. = FALSE)
    } else {
      wb <- openxlsx::createWorkbook()
      excel_row_limit <- 1048576L
      for (nm in names(tables)) {
        openxlsx::addWorksheet(wb, nm)
        tab <- tables[[nm]]
        if (nrow(tab) > excel_row_limit - 1L) {
          tab <- data.frame(
            note = sprintf(paste0("%d rows exceed the Excel sheet limit; ",
                                  "see %s_%s.csv"), nrow(tab), prefix, nm),
            stringsAsFactors = FALSE)
        }
        openxlsx::writeData(wb, nm, tab)
      }
      path <- file.path(dir, sprintf("%s_results.xlsx", prefix))
      openxlsx::saveWorkbook(wb, path, overwrite = overwrite)
      written <- c(written, path)
    }
  }

  invisible(written)
}
