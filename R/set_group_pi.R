#' Create Group Probability of Infection Function
#'
#' Returns a function that computes the average probability of infection
#' between two ages \code{a} and \code{b}, based on the specified model type.
#' This is often used in catalytic models for age-stratified infection data.
#'
#' Supported model types:
#' \itemize{
#'   \item \code{"MuenchGeneral"}: General Muench model with three parameters (\eqn{k}, \eqn{l}, \eqn{\lambda}).
#'   \item \code{"MuenchRestricted"}: Muench model with fixed \eqn{k} and \eqn{l}, estimating only \eqn{\lambda}.
#'   \item \code{"Griffiths"}: Griffiths model with a change point \eqn{\tau} and parameters \eqn{\gamma_0}, \eqn{\gamma_1}.
#'   \item \code{"Farringtons"}: Farrington’s model with parameters \eqn{\gamma_0}, \eqn{\gamma_1}, \eqn{\gamma_2}.
#'   \item \code{"PiecewiseConstant"}: Piecewise constant force of infection model with specified age cutoffs.
#' }
#'
#' @param type Character string specifying the model type. Must be one of:
#'   \code{"MuenchGeneral"}, \code{"MuenchRestricted"}, \code{"Griffiths"},
#'   \code{"Farringtons"}, \code{"PiecewiseConstant"}.
#' @param model_fixed_params List of fixed parameters required by the chosen model type.
#'   For example:
#'   \itemize{
#'     \item MuenchRestricted: \code{list(k = ..., l = ...)}
#'     \item Griffiths: \code{list(tau = ...)}
#'     \item PiecewiseConstant: \code{list(upper_cutoffs = ...)}
#'   }
#'
#' @return A function \code{group_pi(a, b, par)} where:
#'   \itemize{
#'     \item \code{a} is the lower age bound.
#'     \item \code{b} is the upper age bound.
#'     \item \code{par} is a numeric vector of model parameters to be estimated.
#'   }
#'   The returned function computes the average infection probability over \eqn{[a, b]}.
#'
#' @details
#' For some models, the computation involves:
#' \itemize{
#'   \item Analytical formulas (\code{MuenchGeneral}, \code{MuenchRestricted}, \code{Griffiths}).
#'   \item Numerical integration (\code{Farringtons}).
#'   \item Stepwise integration over defined age intervals (\code{PiecewiseConstant}).
#' }
#'
#' The Griffiths model implementation uses the error function via \code{pnorm},
#' and Farrington's model uses \code{\link[stats]{integrate}} with error handling.
#'
#' @importFrom stats pnorm integrate
#' @export
set_group_pi <- function(type, model_fixed_params) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (type == "MuenchGeneral") {
    group_pi <- function(a, b, par) {
      k <- par[1]
      l <- par[2]
      foi <- par[3]
      if (foi < 1e-6) foi <- 1e-6
      return(k*(l-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a))))
    }
  }

  else if (type == "MuenchRestricted") {
    k <- model_fixed_params$k
    l <- model_fixed_params$l
    group_pi <- function(a, b, par) {
      foi <- par[1]
      return(k*(l-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a))))
      # return(k * ((exp(-b * foi) * (b * l * foi * exp(b * foi) + 1)) / foi - (exp(-a * foi) * (a * l * foi * exp(a * foi) + 1)) / foi)) / (b - a)
    }
  }

  else if (type == "Griffiths") {
    tau <- model_fixed_params$tau
    group_pi <- function(a, b, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]

      if (b <= tau) return(0)

      # Adjust bounds if partially below tau
      a_adj <- max(a, tau)
      b_adj <- b

      sqrt_pi <- sqrt(pi)
      sqrt2 <- sqrt(2)
      sqrt_gamma0 <- sqrt(gamma0)

      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

      exp_term <- exp((gamma0 * a_adj^2) / 2 + gamma0 * gamma1 * a_adj + (gamma0 * gamma1^2) / 2)

      erf1 <- erf(sqrt_gamma0 * (sqrt2 * a_adj + sqrt2 * gamma1) / 2)
      erf2 <- erf(sqrt_gamma0 * (sqrt2 * b_adj + sqrt2 * gamma1) / 2)

      A_term <- sqrt_pi * sqrt_gamma0 * sqrt2 * exp_term * (erf1 - erf2)
      B_term <- 2 * gamma0 * (b_adj - a_adj)

      integral <- (A_term + B_term) / (2 * gamma0)

      avg_pi <- integral / (b - a)

      return(avg_pi)
    }
  }

  else if (type == "Farringtons") {
    group_pi <- function(a, b, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
      pi_func <- function(t) {
        return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
      }
      tryCatch(
        integrate(pi_func, a, b)$value / (b - a),
        error = function(e) {
          message(paste("Integration failed on interval [", a, ",", b, "] for gamma0 =",
                        gamma0, ", gamma1 =", gamma1, ", gamma2 =", gamma2, ":", e$message))
          return(NA)
        }
      )
    }
  }

  else if (type == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    group_pi <- function(a, b, par) {
      foi_pieces <- par[ names(par) != "rho" ]
      # warning("Current foi_pieces: ", paste(round(foi_pieces, 4), collapse = ", "))
      interval_lengths_a <- pmax(pmin(a, upper_cutoffs) - lower_cutoffs, 0)
      cum_foi_a <- sum(foi_pieces * interval_lengths_a)
      pi_a <- 1-exp(-cum_foi_a)
      interval_lengths_b <- pmax(pmin(b, upper_cutoffs) - lower_cutoffs, 0)
      cum_foi_b <- sum(foi_pieces * interval_lengths_b)
      pi_b <- 1-exp(-cum_foi_b)
      return((pi_b-pi_a)/(b-a))
    }
  }


  else {
    return(NA)
  }

  return(group_pi)
}
