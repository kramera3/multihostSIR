#' Required columns for a user-supplied community
#'
#' A complete community table has one row per species and these columns. When a
#' user supplies all of them, together with a transmission matrix, the model is
#' fully determined and needs no random draws.
#' @noRd
required_community_columns <- function() {
  c("type", "mass", "mu", "gamma", "alpha", "N", "competence")
}

#' Read a host community from a data frame, CSV, or Excel file
#'
#' Turns a species table into the community object the rest of the package uses.
#' The input may be a data frame already in memory, a path to a `.csv` file, or a
#' path to an Excel `.xlsx`/`.xls` file.
#'
#' @details
#' One row per species. The columns are:
#'
#' \tabular{lll}{
#'   **Column** \tab **Meaning** \tab **Needed?** \cr
#'   `type`       \tab `"Prey"` or `"Predator"`          \tab always \cr
#'   `mass`       \tab body mass (kg)                    \tab always \cr
#'   `mu`         \tab background mortality (per day)    \tab for a complete table \cr
#'   `gamma`      \tab recovery rate (per day)           \tab for a complete table \cr
#'   `alpha`      \tab disease-induced mortality (per day) \tab for a complete table \cr
#'   `N`          \tab density (individuals per km^2)    \tab for a complete table \cr
#'   `competence` \tab within-species transmission \eqn{\beta} \tab for a complete table \cr
#'   `species`    \tab a name or label                  \tab optional \cr
#' }
#'
#' If a rate column is absent it is filled from body mass using `allo`, exactly
#' as [build_community()] would. A table is "complete" only when every column
#' above is present and supplied for every species. Reading Excel needs the
#' \pkg{readxl} package.
#'
#' @param x A data frame, or a path to a `.csv`, `.xlsx`, or `.xls` file.
#' @param allo An allometry object from [allometry()], used to fill any missing
#'   rate columns.
#' @param sheet For Excel input, the sheet name or number. Defaults to the first
#'   sheet.
#'
#' @return A community data frame with the columns [build_community()] produces,
#'   carrying an attribute `"complete"` that records whether every rate was
#'   supplied by the user (`TRUE`) or any were filled from allometry (`FALSE`).
#'
#' @examples
#' # From a data frame already in memory
#' species <- data.frame(
#'   species    = c("mouse", "deer", "fox"),
#'   type       = c("Prey", "Prey", "Predator"),
#'   mass       = c(0.02, 50, 5),
#'   mu         = c(0.013, 0.002, 0.003),
#'   gamma      = c(0.27, 0.04, 0.07),
#'   alpha      = c(0.027, 0.004, 0.007),
#'   N          = c(9000, 27, 6),
#'   competence = c(0.5, 0.1, 0.2))
#' comm <- read_community(species)
#' attr(comm, "complete")   # TRUE: nothing was filled in
#'
#' \dontrun{
#' # From files
#' comm <- read_community("my_species.csv")
#' comm <- read_community("my_species.xlsx", sheet = "communities")
#' }
#'
#' @seealso [read_transmission_matrix()], [run_fixed_community()],
#'   [build_community()]
#' @export
read_community <- function(x, allo = allometry(), sheet = NULL) {
  if (!inherits(allo, "multihost_allometry")) {
    stop("`allo` must be created by allometry().", call. = FALSE)
  }

  df <- read_table_input(x, sheet = sheet, what = "community")

  if (!"type" %in% names(df)) {
    stop("The community table needs a `type` column (\"Prey\"/\"Predator\").",
         call. = FALSE)
  }
  if (!"mass" %in% names(df) &&
      !all(c("mu", "gamma", "alpha", "N") %in% names(df))) {
    stop(paste0("The community table needs a `mass` column, or all of `mu`, ",
                "`gamma`, `alpha` and `N` so that mass is not required."),
         call. = FALSE)
  }
  df$type <- as.character(df$type)
  if (!all(df$type %in% c("Prey", "Predator"))) {
    bad <- unique(df$type[!df$type %in% c("Prey", "Predator")])
    stop(sprintf("`type` must be \"Prey\" or \"Predator\"; found: %s.",
                 paste(bad, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(df)

  # Fill any missing rate columns from body mass, and record whether we had to.
  filled <- character(0)
  if ("mass" %in% names(df)) {
    scaled <- apply_allometry(df$mass, allo)
    for (col in c("mu", "gamma", "alpha", "N")) {
      if (!col %in% names(df) || anyNA(df[[col]])) {
        df[[col]] <- scaled[[col]]
        filled <- c(filled, col)
      }
    }
  }

  if (!"competence" %in% names(df) || anyNA(df$competence)) {
    # Fall back to the midpoint of the default medium-competence range, so a
    # single deterministic value is used rather than a random draw.
    df$competence <- 0.2
    filled <- c(filled, "competence")
  }

  if (is.null(df$species)) df$species <- seq_len(n)
  df$trophic_weight <- ifelse(df$type == "Predator", 2, 1)

  # competence_category is only used for reporting; derive a label from the
  # supplied competence so downstream summaries still work.
  df$competence_category <- factor(
    ifelse(df$competence < 0.1, "Low",
           ifelse(df$competence < 0.3, "Med", "High")),
    levels = c("Low", "Med", "High"))

  df$beta_ii <- df$competence / df$N

  keep <- c("species", "type", "mass", "trophic_weight",
            "mu", "gamma", "alpha", "N",
            "competence_category", "competence", "beta_ii")
  df <- df[, keep, drop = FALSE]

  attr(df, "complete") <- length(filled) == 0L
  attr(df, "filled") <- filled
  df
}

#' Read a transmission matrix from a matrix, CSV, or Excel file
#'
#' Reads an \eqn{n \times n} matrix of between-species transmission rates. Row
#' \eqn{i}, column \eqn{j} is the rate at which one infected individual of
#' species \eqn{j} infects susceptibles of species \eqn{i}.
#'
#' @param x A numeric matrix, or a path to a `.csv`, `.xlsx`, or `.xls` file
#'   holding a square grid of numbers. A header row and first-column labels are
#'   allowed and are used to check the ordering against the community.
#' @param community Optional community (from [read_community()] or
#'   [build_community()]). When given, the matrix dimension is checked against
#'   the number of species, and any row/column labels are checked against the
#'   community's species names.
#' @param sheet For Excel input, the sheet name or number.
#'
#' @return A numeric \eqn{n \times n} matrix.
#'
#' @examples
#' m <- matrix(c(0.05, 0.01, 0.00,
#'               0.01, 0.03, 0.02,
#'               0.00, 0.02, 0.04), nrow = 3, byrow = TRUE)
#' beta <- read_transmission_matrix(m)
#' beta
#'
#' \dontrun{
#' beta <- read_transmission_matrix("my_matrix.xlsx")
#' }
#'
#' @seealso [read_community()], [run_fixed_community()], [transmission_matrix()]
#' @export
read_transmission_matrix <- function(x, community = NULL, sheet = NULL) {
  if (is.matrix(x)) {
    m <- x
  } else {
    df <- read_table_input(x, sheet = sheet, what = "matrix",
                           allow_rownames = TRUE)
    # Drop a leading label column if the first column is non-numeric.
    if (!is.numeric(df[[1L]])) {
      rn <- as.character(df[[1L]])
      df <- df[, -1L, drop = FALSE]
      m <- as.matrix(df)
      rownames(m) <- rn
    } else {
      m <- as.matrix(df)
    }
  }

  storage.mode(m) <- "double"
  if (nrow(m) != ncol(m)) {
    stop(sprintf("The transmission matrix must be square; got %d x %d.",
                 nrow(m), ncol(m)), call. = FALSE)
  }
  if (anyNA(m) || any(!is.finite(m)) || any(m < 0)) {
    stop("Transmission-matrix entries must be finite and non-negative.",
         call. = FALSE)
  }

  if (!is.null(community)) {
    check_community(community)
    if (nrow(m) != nrow(community)) {
      stop(sprintf(paste0("The transmission matrix is %d x %d but the ",
                          "community has %d species."),
                   nrow(m), ncol(m), nrow(community)), call. = FALSE)
    }
    labels <- rownames(m) %||% colnames(m)
    sp <- as.character(community$species %||% seq_len(nrow(community)))
    if (!is.null(labels) && !identical(as.character(labels), sp) &&
        !setequal(as.character(labels), sp)) {
      warning(paste0("Matrix row/column labels do not match the community's ",
                     "species names; assuming they are in the same order."),
              call. = FALSE)
    }
    dimnames(m) <- list(sp, sp)
  }
  m
}

#' Run one fixed, fully specified community (no iteration)
#'
#' When a community and its transmission matrix are both supplied, the epidemic
#' is deterministic: there is nothing random to average over, so the model is
#' solved exactly once. This is the counterpart to [run_multihost_sir()], which
#' repeats random draws.
#'
#' @param community A community from [read_community()] or [build_community()].
#' @param beta A transmission matrix from [read_transmission_matrix()], or any
#'   \eqn{n \times n} numeric matrix. If `NULL`, the matrix is built from the
#'   community with [transmission_matrix()] and `coefficients` (in which case the
#'   run is still deterministic, but uses the package's transmission rules rather
#'   than your own numbers).
#' @param duration Simulation horizon in days.
#' @param coefficients Trophic coefficients, used only when `beta = NULL`.
#' @param ... Passed to [simulate_epidemic()].
#'
#' @return An object of class `"multihost_fixed"`: a list with `results` (a
#'   one-row data frame of community-level outcomes), `species` (per-species
#'   metrics), `trajectory` (daily prevalence), `R0`, `community`, `beta`, and
#'   `assumptions` (a record of what was supplied versus filled in).
#'
#' @examples
#' species <- data.frame(
#'   type       = c("Prey", "Prey", "Predator"),
#'   mass       = c(0.02, 50, 5),
#'   mu         = c(0.013, 0.002, 0.003),
#'   gamma      = c(0.27, 0.04, 0.07),
#'   alpha      = c(0.027, 0.004, 0.007),
#'   N          = c(9000, 27, 6),
#'   competence = c(0.5, 0.1, 0.2))
#' comm <- read_community(species)
#' beta <- read_transmission_matrix(matrix(c(
#'   0.05, 0.01, 0.00,
#'   0.01, 0.03, 0.00,
#'   0.02, 0.02, 0.04), nrow = 3, byrow = TRUE), community = comm)
#'
#' fit <- run_fixed_community(comm, beta, duration = 200)
#' fit
#' fit$results[, c("richness", "R0", "peak_prevalence")]
#'
#' @seealso [run_multihost_sir()], [read_community()],
#'   [read_transmission_matrix()]
#' @export
run_fixed_community <- function(community, beta = NULL, duration = 365,
                                coefficients = trophic_coefficients(), ...) {
  check_community(community)
  check_scalar_num(duration, "duration", lo = 1, integer = TRUE)

  supplied_beta <- !is.null(beta)
  if (is.null(beta)) {
    beta <- transmission_matrix(community, coefficients = coefficients)
  } else if (!is.matrix(beta) || nrow(beta) != nrow(community) ||
             ncol(beta) != nrow(community)) {
    stop(sprintf("`beta` must be a %d x %d matrix for this community.",
                 nrow(community), nrow(community)), call. = FALSE)
  }

  sim <- simulate_epidemic(community, beta = beta, duration = duration, ...)
  cm <- community_metrics(sim)

  # Keep the per-species infected curves as well, so plots can show each
  # species, not only the community total.
  n <- nrow(community)
  I_cols <- paste0("I", seq_len(n))
  species_curves <- data.frame(
    time = sim$out$time,
    sim$out[, I_cols, drop = FALSE],
    check.names = FALSE)
  sp_names <- as.character(community$species %||% seq_len(n))
  names(species_curves) <- c("time", sp_names)

  complete <- isTRUE(attr(community, "complete")) && supplied_beta
  filled <- attr(community, "filled") %||% character(0)

  assumptions <- list(
    deterministic   = complete,
    beta_supplied   = supplied_beta,
    filled_columns  = filled)

  structure(
    list(results        = cm$row,
         species        = species_metrics(sim),
         trajectory     = data.frame(time = cm$time,
                                     prevalence = cm$prevalence),
         species_curves = species_curves,
         R0             = sim$R0,
         community      = community,
         beta           = beta,
         duration       = as.integer(duration),
         assumptions    = assumptions),
    class = "multihost_fixed")
}

#' @param x A `"multihost_fixed"` object.
#' @param ... Ignored.
#' @rdname run_fixed_community
#' @export
print.multihost_fixed <- function(x, ...) {
  cat("<multihost_fixed: one fully specified community>\n")
  cat(sprintf("  species        : %d\n", nrow(x$community)))
  cat(sprintf("  horizon        : %d days\n", x$duration))
  cat(sprintf("  community R0    : %.3f\n", x$R0))
  cat(sprintf("  peak prevalence : %.1f infected individuals\n",
              x$results$peak_prevalence))
  cat(sprintf("  day of peak     : %g\n", x$results$days_to_peak))
  cat("\n")
  print_assumptions(x$assumptions, x)
  cat("\nUse $results, $species and $trajectory for the numbers.\n")
  invisible(x)
}

#' Read a data frame from memory, CSV, or Excel
#' @noRd
read_table_input <- function(x, sheet = NULL, what = "table",
                             allow_rownames = FALSE) {
  if (is.data.frame(x)) return(x)
  if (!is.character(x) || length(x) != 1L) {
    stop(sprintf("`%s` must be a data frame or a single file path.", what),
         call. = FALSE)
  }
  if (!file.exists(x)) {
    stop(sprintf("File not found: %s", x), call. = FALSE)
  }
  ext <- tolower(tools::file_ext(x))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop(paste0("Reading Excel files needs the 'readxl' package. Install it ",
                  "with install.packages(\"readxl\"), or save your file as CSV ",
                  "and read that instead."), call. = FALSE)
    }
    sheet <- sheet %||% 1L
    as.data.frame(readxl::read_excel(x, sheet = sheet),
                  stringsAsFactors = FALSE)
  } else if (ext == "csv") {
    utils::read.csv(x, stringsAsFactors = FALSE,
                    check.names = FALSE)
  } else {
    stop(sprintf("Unsupported file type '.%s'; use .csv, .xlsx or .xls.", ext),
         call. = FALSE)
  }
}
