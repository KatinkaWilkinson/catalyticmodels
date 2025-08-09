#' Create Group-Averaged Force of Infection Function
#'
#' Returns a function \code{group_foi(a, b, par)} (or \code{group_foi(a, b, spline_pi_t)} for splines)
#' that computes the **average force of infection** over an age interval \eqn{[a, b)}, based on the
#' specified model type.
#'
#' Supported model types:
#' \itemize{
#'   \item \code{"MuenchGeneral"} – constant FOI from the third parameter.
#'   \item \code{"MuenchRestricted"} – constant FOI from the first parameter.
#'   \item \code{"Griffiths"} – FOI is zero until \eqn{\tau}, then increases linearly:
#'         \eqn{\gamma_0 (t + \gamma_1)}; integrates exactly over \eqn{[a,b)}.
#'   \item \code{"Farringtons"} – exponential decay plus constant term; integrates exactly over \eqn{[a,b)}.
#'   \item \code{"PiecewiseConstant"} – constant within age intervals defined by \code{upper_cutoffs};
#'         computes a weighted average for overlaps with \eqn{[a,b)}.
#'   \item \code{"Splines"} – computes average FOI by integrating
#'         \eqn{\lambda(t) = \pi'(t)/(1 - \pi(t))} from a smoothed prevalence spline, then dividing by interval width.
#' }
#'
#' @param type Character string specifying the model type.
#' @param model_fixed_params Optional named list of fixed parameters required by certain models:
#'   \itemize{
#'     \item \code{"Griffiths"}: \code{list(tau = ...)}
#'     \item \code{"PiecewiseConstant"}: \code{list(upper_cutoffs = ...)}
#'   }
#'
#' @return A function \code{group_foi(a, b, par)} (or \code{group_foi(a, b, spline_pi_t)} for splines), where:
#' \itemize{
#'   \item \code{a} is the lower bound of the age interval.
#'   \item \code{b} is the upper bound of the age interval.
#'   \item \code{par} is a numeric vector of model parameters (for spline models, pass a spline object instead).
#'   \item The returned value is the **average** FOI over \eqn{[a, b)}.
#' }
#'
#' @details
#' \strong{Piecewise models:} For \code{"PiecewiseConstant"}, overlaps between each defined interval and \eqn{[a,b)} are
#' computed, and a weighted average is returned.
#'
#' \strong{Spline models:} The function uses numerical integration via \code{\link[stats]{integrate}} to approximate
#' the mean FOI, with safeguards against division by zero.
#'
#' Base R functions like \code{ifelse}, \code{pmin}, \code{pmax}, \code{findInterval}, and \code{integrate} are used,
#' so no package imports are required.
#'
#' @export
set_group_foi <- function(type, model_fixed_params = NA) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (type == "MuenchGeneral") {
    group_foi <- function(a, b, par) {
      foi <- par[3]
      return(foi) # constant foi
    }
  }

  else if (type == "MuenchRestricted") {
    group_foi <- function(a, b, par) {
      foi <- par[1]
      return(foi) # constant foi
    }
  }

  else if (type == "Griffiths") {
    tau <- model_fixed_params$tau
    group_foi <- function(a, b, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]

      if (b <= tau) return(0)

      # Adjust bounds if partially below tau
      a_adj <- max(a, tau)
      b_adj <- b
      integral <- ((b_adj - a_adj) * (2 * gamma0 * gamma1 + (b_adj + a_adj) * gamma0)) / 2
      avg_foi <- integral / (b - a)
      return(avg_foi)
    }
  }

  else if (type == "Farringtons") {
    group_foi <- function(a, b, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
      integral <- (((-gamma1 - b * gamma0) * gamma2 - gamma0) * exp(-b * gamma2) + exp(-a * gamma2) * ((b - a) * gamma1 * gamma2^2 * exp(a * gamma2) + (gamma1 + a * gamma0) * gamma2 + gamma0)) / gamma2^2
      avg_foi <- integral / (b - a)
      return(avg_foi)
    }
  }

  else if (type == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    group_foi <- function(a, b, par) {
      # Determine overlap between each piecewise interval and [a, b)
      overlaps <- pmin(b, upper_cutoffs) - pmax(a, lower_cutoffs)
      overlaps[overlaps < 0] <- 0  # No overlap = set to zero
      total_length <- b - a
      weighted_avg <- sum(par * overlaps) / total_length
      return(weighted_avg)
    }
  }

  else if (type == "Splines") {
    group_foi <- function(a, b, spline_pi_t) {
      integrand <- function(x) {
        pi_x <- predict(spline_pi_t, x)$y
        dpi_x <- predict(spline_pi_t, x, deriv = 1)$y
        pi_x <- pmin(pi_x, 1 - 1e-8) # to prevent a divide by 0 error
        dpi_x / (1 - pi_x) # estimated foi
      }
      result <- tryCatch({
        integrate(integrand, a, b)$value / (b - a) #  approximate the mean FOI over [a, b)  by integrating the FOI (π′ / (1−π)) over that interval and dividing by the interval width
      }, error = function(e) NA)
      return(result)
    }
  }

  else {
    return(NA)
  }

  return(group_foi)
}
