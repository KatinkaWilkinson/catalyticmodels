#' Create Instantaneous Force of Infection Function
#'
#' Returns a function \code{foi_t(t, par)} that computes the instantaneous force
#' of infection at age \code{t}, based on the specified model type.
#'
#' Supported model types:
#' \itemize{
#'   \item \code{"MuenchGeneral"} – constant FOI from the third parameter.
#'   \item \code{"MuenchRestricted"} – constant FOI from the first parameter.
#'   \item \code{"Griffiths"} – FOI is zero until \eqn{\tau}, then increases linearly: \eqn{\gamma_0 (t + \gamma_1)}.
#'   \item \code{"Farringtons"} – decaying exponential plus constant term.
#'   \item \code{"PiecewiseConstant"} – constant within age intervals defined by \code{upper_cutoffs}.
#'   \item \code{"Splines"} – computed from a smoothed prevalence spline, using:
#'     \eqn{\lambda(t) = \pi'(t) / \left( 1 - \pi(t) \right)}.
#' }
#'
#' @param type Character string specifying the model type.
#' @param model_fixed_params Optional named list of fixed parameters required by certain models:
#'   \itemize{
#'     \item \code{"Griffiths"}: \code{list(tau = ...)}
#'     \item \code{"PiecewiseConstant"}: \code{list(upper_cutoffs = ...)}
#'   }
#'
#' @return A function \code{foi_t(t, par)} (or \code{foi_t(t, spline_pi_t)} for splines), where:
#' \itemize{
#'   \item \code{t} is age (numeric vector).
#'   \item \code{par} is a numeric vector of model parameters.
#'   \item The returned value is the instantaneous FOI at each age in \code{t}.
#' }
#' If \code{foi_functional_form} is unrecognised, returns \code{NA}.
#'
#' @details
#' For spline-based models, FOI is calculated as the derivative of the prevalence spline divided by
#' \eqn{1 - \pi(t)} (capped to avoid division by zero). This ensures the FOI is well-defined for all \code{t}.
#' Base R functions like \code{ifelse}, \code{findInterval}, and \code{pmin} are used, so no special imports are required.
#'
#' @export
set_foi_t <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA) { # foi_functional_form is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (foi_functional_form == "Constant") {
    foi_t <- function(t, par) {
      foi <- par[["foi"]]
      return(foi) # constant foi.
    }
  }

  else if (foi_functional_form == "Griffiths") {
    tau <- model_fixed_params$tau
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
      return((gamma0*t + gamma1)*exp(-gamma2*t)+gamma1)
    }
  }

  else if (foi_functional_form == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    foi_t <- function(t, par) {
      # Each interval is [lower_cutoffs[i], upper_cutoffs[i])
      interval_index <- findInterval(t, lower_cutoffs, rightmost.closed = FALSE)
      return(par[interval_index])
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

