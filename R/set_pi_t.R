#' Create Prevalence Function \eqn{\pi(t)} for a Given Catalytic Model Type
#'
#' Returns a function \eqn{\pi(t, \theta)} that computes the modelled prevalence
#' at age \code{t} for the specified catalytic model type and parameters.
#'
#' @param type Character string specifying the model type. Supported types are:
#' \itemize{
#'   \item \code{"MuenchGeneral"} – General Muench model with \eqn{k}, \eqn{l}, and FOI as parameters.
#'   \item \code{"MuenchRestricted"} – Muench model with fixed \eqn{k} and \eqn{l}, estimating FOI only.
#'   \item \code{"Griffiths"} – Griffiths model with FOI defined piecewise before and after \eqn{\tau}.
#'   \item \code{"Farringtons"} – Farrington’s model with parameters \eqn{\gamma_0}, \eqn{\gamma_1}, \eqn{\gamma_2}.
#'   \item \code{"PiecewiseConstant"} – Piecewise-constant FOI model defined over fixed age intervals.
#' }
#' @param model_fixed_params Optional named list of fixed parameters required by some models:
#' \itemize{
#'   \item For \code{"MuenchRestricted"}: \code{k}, \code{l}.
#'   \item For \code{"Griffiths"}: \code{tau}.
#'   \item For \code{"PiecewiseConstant"}: \code{upper_cutoffs} (numeric vector of upper age bounds).
#' }
#'
#' @return A function \code{pi_t(t, par)} where:
#' \itemize{
#'   \item \code{t} is a numeric vector of ages.
#'   \item \code{par} is a numeric vector of model parameters to be estimated.
#'   \item The returned value is the modelled prevalence at each age in \code{t}.
#' }
#' If the \code{type} is not recognised, returns \code{NA}.
#'
#' @details
#' The returned function implements the prevalence model \eqn{\pi(t)} corresponding to the
#' specified catalytic model type. For \code{"PiecewiseConstant"}, prevalence is computed
#' as \eqn{1 - \exp(-\sum \lambda_i \cdot L_i)}, where \eqn{\lambda_i} is the FOI in each interval
#' and \eqn{L_i} is the overlap of age \code{t} with that interval.
#'
#' @export
set_pi_t <- function(type, model_fixed_params = NA) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (type == "MuenchGeneral") {
    pi_t <- function(t, par) {
      k <- par[1]
      l <- par[2]
      foi <- par[3]
      return(k * (l - exp(-foi * t)))
    }
  }

  else if (type == "MuenchRestricted") {
    k <- model_fixed_params$k
    l <- model_fixed_params$l
    pi_t <- function(t, par) {
      foi <- par[1]
      return(k * (l - exp(-foi * t)))
    }
  }

  else if (type == "Griffiths") {
    tau <- model_fixed_params$tau
    pi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      return(ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau)))))
    }
  }

  else if (type == "Farringtons") {
    pi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
    }
  }

  else if (type == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])

    pi_t <- function(t, par) {
      foi_pieces <- par[ names(par) != "rho" ]
      interval_lengths <- pmax(pmin(t, upper_cutoffs) - lower_cutoffs, 0)
      cum_foi <- sum(foi_pieces * interval_lengths)
      return(1 - exp(-cum_foi))
    }
  }

  else {
    return(NA)
  }

  return(pi_t)
}
