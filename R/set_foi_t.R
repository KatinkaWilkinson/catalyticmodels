#' Create an FOI function for a chosen model (helper)
#'
#' Returns a function that computes the instantaneous force of infection (FOI)
#' at age \code{t} for the specified \code{foi_functional_form}. This is a
#' helper used internally by \code{\link{FoiFromCatalyticModel}} and related
#' wrappers; it is not intended for direct end-user calls.
#'
#' Supported forms:
#' \itemize{
#'   \item \code{"Constant"} — FOI is constant (uses \code{par[["foi"]]}), set to 0 for ages < \code{tau}.
#'   \item \code{"Linear"} — FOI is \code{m * t + c} (uses \code{par[["m"]]}, \code{par[["c"]]}), 0 for ages < \code{tau}.
#'   \item \code{"Griffiths"} — FOI is 0 up to \eqn{\tau} then \eqn{\gamma_0 (t + \gamma_1)}; requires \code{model_fixed_params$tau}.
#'   \item \code{"Farringtons"} — \eqn{( \gamma_0 t - \gamma_1 ) e^{-\gamma_2 t} + \gamma_1}, 0 for ages < \code{tau}.
#'   \item \code{"PiecewiseConstant"} — Constant within intervals \[a, b), from \code{model_fixed_params$upper_cutoffs}.
#'   \item \code{"Splines"} + \code{catalytic_model_type == "SimpleCatalytic"} — \eqn{\lambda(t) = \pi'(t)/(1-\pi(t))}.
#'   \item \code{"Splines"} + \code{catalytic_model_type == "WaningImmunity"} — \eqn{\lambda(t) = \{\pi'(t)+w\,\pi(t)\}/(1-\pi(t))}, requires \code{model_fixed_params$w}.
#' }
#'
#' @param catalytic_model_type Character. Model family (e.g., \code{"SimpleCatalytic"}, \code{"WaningImmunity"})
#'   used to disambiguate spline variants.
#' @param foi_functional_form Character. One of \code{"Constant"}, \code{"Linear"}, \code{"Griffiths"},
#'   \code{"Farringtons"}, \code{"PiecewiseConstant"}, \code{"Splines"}.
#' @param model_fixed_params Optional named list for required settings in some forms:
#'   \itemize{
#'     \item \code{Griffiths}: \code{list(tau = ...)}
#'     \item \code{PiecewiseConstant}: \code{list(upper_cutoffs = numeric())} — ascending upper bounds per interval
#'     \item \code{WaningImmunity + Splines}: \code{list(w = ...)}
#'   }
#' @param tau Numeric. Age threshold; FOI is set to 0 for ages < \code{tau} in \code{Constant}, \code{Linear},
#'   \code{Farringtons}, and \code{PiecewiseConstant}. For \code{Griffiths}, \code{model_fixed_params$tau} is used.
#'
#' @return A function:
#' \itemize{
#'   \item Analytic forms: \code{foi_t(t, par)} where \code{par} is a named vector.
#'   \item Spline forms: \code{foi_t(t, spline_pi_t)} where \code{spline_pi_t} is a \code{\link[stats]{smooth.spline}} fit.
#' }
#' If the combination is unrecognised, returns \code{NA}.
#'
#' @details
#' Spline forms use \code{\link[stats]{predict}} for \eqn{\pi(t)} and \eqn{\pi'(t)}; \eqn{\pi(t)} is capped below 1 to
#' avoid division by zero. For \code{PiecewiseConstant}, intervals are treated as \code{[a, b)} and \code{par} must
#' supply one FOI per interval in the same order as \code{upper_cutoffs}.
#'
#' @note Helper for \code{\link{FoiFromCatalyticModel}}; not intended for independent use.
#'       Consider not exporting this function. If you keep it exported for advanced users,
#'       treat it as internal API subject to change.
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}, \code{\link{FoiFromCatalyticModel_unparallelised}}
#'
#' @keywords internal
#'
#' @importFrom stats predict
#' @noRd
set_foi_t <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, tau = 0) { # foi_functional_form is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (foi_functional_form == "Constant") {
    foi_t <- function(t, par) {
      foi <- par[["foi"]]
      return(ifelse(t<tau, 0, foi)) # constant foi
    }
  }

  else if (foi_functional_form == "Linear") {
    foi_t <- function(t, par) {
      m <- par[["m"]]
      c <- par[["c"]]
      return(ifelse(t<tau, 0, m*t+c))
    }
  }

  else if (foi_functional_form == "Griffiths") {
    tau <- model_fixed_params$tau
    force(tau)
    foi_t <- function(t, par) {
      gamma0 <- par[["gamma0"]]
      gamma1 <- par[["gamma1"]]
      return(ifelse(t <= tau, 0, gamma0 * (t + gamma1)))
    }
    return(foi_t)
  }

  else if (foi_functional_form == "Farringtons") {
    foi_t <- function(t, par) {
      gamma0 <- par[["gamma0"]]
      gamma1 <- par[["gamma1"]]
      gamma2 <- par[["gamma2"]]
      return(ifelse(t<tau, 0, (gamma0*t - gamma1)*exp(-gamma2*t)+gamma1))
    }
  }

  else if (foi_functional_form == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    force(upper_cutoffs)
    force(lower_cutoffs)
    foi_t <- function(t, par) {
      # Each interval is [lower_cutoffs[i], upper_cutoffs[i])
      interval_index <- findInterval(t, lower_cutoffs, rightmost.closed = FALSE)
      return(ifelse(t<tau, 0, par[interval_index]))
    }
  }

  else if (!is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Splines") {
    foi_t <- function(t, spline_pi_t) {
      pi <- predict(spline_pi_t, t)$y # y is the smooth.spline response name. "y" here is predicting pi(t) (we inputted y/n - this is what we are trying to predict. Remember, y/n=pi(t))
      dpi_t_dt <- predict(spline_pi_t, t, deriv=1)$y
      pi <- pmin(pi, 1 - 1e-8) # to prevent a divide by 0 error
      foi <- dpi_t_dt / (1 - pi) # dS_dt = -foi*(1-pi), therefore foi = -dS_dt/(1-pi) = dpi_dt/(1-pi)

      # dpi_t_dt <- diff(pi_t)/diff(t)
      # foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)]) this is what I used to have but apparently using splines is better for finding the derivative of pi_t

      return(foi)
    }
  }

  else if (!is.na(catalytic_model_type) && catalytic_model_type == "WaningImmunity" && foi_functional_form == "Splines") {
    w <- model_fixed_params$w
    force(w)
    foi_t <- function(t, spline_pi_t) {
      pi <- predict(spline_pi_t, t)$y # y is the smooth.spline response name. "y" here is predicting pi(t) (we inputted y/n - this is what we are trying to predict. Remember, y/n=pi(t))
      dpi_t_dt <- predict(spline_pi_t, t, deriv=1)$y
      pi <- pmin(pi, 1 - 1e-8) # to prevent a divide by 0 error
      foi <- (dpi_t_dt + pi*w) / (1 - pi)

      # dpi_t_dt <- diff(pi_t)/diff(t)
      # foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)]) this is what I used to have but apparently using splines is better for finding the derivative of pi_t

      return(foi)
    }
  }

  else {
    return(NA)
  }

  return(foi_t)
}

