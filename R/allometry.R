#' Allometric scaling coefficients
#'
#' Collects the intercepts and exponents of the body-mass power laws
#' \eqn{Y = a M^{b}} used to derive demographic and epidemiological rates from
#' body mass \eqn{M} (kg). The defaults reproduce the quarter-power scaling used
#' throughout the package.
#'
#' @details
#' The four scaled quantities are:
#'
#' \tabular{lll}{
#'   **Symbol** \tab **Meaning** \tab **Default** \cr
#'   \eqn{\mu} \tab background mortality rate (day\eqn{^{-1}}) \tab
#'     \eqn{0.005 M^{-0.25}} \cr
#'   \eqn{\gamma} \tab recovery rate (day\eqn{^{-1}}) \tab
#'     \eqn{0.1 M^{-0.25}} \cr
#'   \eqn{N} \tab carrying-capacity density (ind km\eqn{^{-2}}) \tab
#'     \eqn{500 M^{-0.75}} \cr
#'   \eqn{\alpha} \tab disease-induced mortality (day\eqn{^{-1}}) \tab
#'     \eqn{0.01 M^{-0.25}}
#' }
#'
#' Biological justification for the exponents:
#'
#' * Mortality scales as \eqn{M^{-0.25}}: quarter-power scaling of metabolic and
#'   physiological rates, and of lifespan (Kleiber 1932; Peters 1983;
#'   Speakman 2005; Healy et al. 2014). Large animals live longer.
#' * Recovery scales as \eqn{M^{-0.25}}: physiological process *durations* scale
#'   as \eqn{M^{0.25}}, so the corresponding *rates* scale as \eqn{M^{-0.25}}
#'   (Lindstedt and Calder 1981; Peters 1983). Large animals recover more slowly.
#' * Density scales as \eqn{M^{-0.75}} (Damuth's law; Damuth 1981, 1987). Large
#'   animals are rarer.
#'
#' The unit of time throughout the package is one day.
#'
#' @param a_mort,b_mort Intercept and exponent for background mortality.
#' @param a_recovery,b_recovery Intercept and exponent for the recovery rate.
#' @param a_density,b_density Intercept and exponent for carrying-capacity
#'   density.
#' @param a_virulence,b_virulence Intercept and exponent for disease-induced
#'   mortality.
#'
#' @return An object of class `"multihost_allometry"`: a list holding the eight
#'   coefficients.
#'
#' @examples
#' allo <- allometry()
#' allo
#'
#' # A 0.02 kg host has a much higher per-capita mortality than a 50 kg host
#' allo$a_mort * c(0.02, 50)^allo$b_mort
#'
#' # Slower recovery, doubled virulence
#' allometry(a_recovery = 0.05, a_virulence = 0.02)
#'
#' @seealso [build_community()], [plot_allometry()]
#' @export
allometry <- function(a_mort      = 0.005, b_mort      = -0.25,
                      a_recovery  = 0.1,   b_recovery  = -0.25,
                      a_density   = 500,   b_density   = -0.75,
                      a_virulence = 0.01,  b_virulence = -0.25) {

  check_scalar_num(a_mort,      "a_mort",      lo = 0)
  check_scalar_num(a_recovery,  "a_recovery",  lo = 0)
  check_scalar_num(a_density,   "a_density",   lo = 0)
  check_scalar_num(a_virulence, "a_virulence", lo = 0)
  check_scalar_num(b_mort,      "b_mort")
  check_scalar_num(b_recovery,  "b_recovery")
  check_scalar_num(b_density,   "b_density")
  check_scalar_num(b_virulence, "b_virulence")

  structure(
    list(a_mort      = a_mort,      b_mort      = b_mort,
         a_recovery  = a_recovery,  b_recovery  = b_recovery,
         a_density   = a_density,   b_density   = b_density,
         a_virulence = a_virulence, b_virulence = b_virulence),
    class = "multihost_allometry")
}

#' @param x An object of class `"multihost_allometry"`.
#' @param ... Ignored.
#' @rdname allometry
#' @export
print.multihost_allometry <- function(x, ...) {
  cat("<multihost allometry>\n")
  rows <- data.frame(
    parameter = c("mortality (mu)", "recovery (gamma)",
                  "density (N)", "virulence (alpha)"),
    intercept = c(x$a_mort, x$a_recovery, x$a_density, x$a_virulence),
    exponent  = c(x$b_mort, x$b_recovery, x$b_density, x$b_virulence),
    stringsAsFactors = FALSE)
  print(rows, row.names = FALSE)
  cat("Rates are per day; density is individuals per square kilometre.\n")
  invisible(x)
}

#' Apply allometric scaling to a vector of body masses
#'
#' @param mass Numeric vector of body masses (kg).
#' @param allo An object from [allometry()].
#' @return A data frame with columns `mu`, `gamma`, `N` and `alpha`.
#' @noRd
apply_allometry <- function(mass, allo) {
  data.frame(
    mu    = allo$a_mort      * mass^allo$b_mort,
    gamma = allo$a_recovery  * mass^allo$b_recovery,
    N     = allo$a_density   * mass^allo$b_density,
    alpha = allo$a_virulence * mass^allo$b_virulence)
}
