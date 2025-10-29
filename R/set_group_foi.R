#' Create a group-FOI (mean over interval) function (helper)
#'
#' Returns a function that computes the mean force of infection (FOI) over an
#' age interval \[a, b), i.e. \eqn{\bar{\lambda}_{[a,b)}}, for the specified
#' \code{foi_functional_form}. This is a helper used internally by
#' \code{\link{FoiFromCatalyticModel}} and related wrappers; it is not intended
#' for direct end-user calls.
#'
#' Supported forms:
#' \itemize{
#'   \item \code{"Constant"} — mean FOI equals \code{par[["foi"]]}.
#'   \item \code{"Linear"} — \code{m * t + c}; mean over \[a, b) is the average of the endpoints:
#'         \code{(m*a + c + m*b + c)/2}.
#'   \item \code{"Griffiths"} — 0 if \code{b <= tau}; otherwise the mean of
#'         \eqn{\gamma_0 (t + \gamma_1)} over \[a, b). Requires \code{model_fixed_params$tau}.
#'   \item \code{"Farringtons"} — mean of
#'         \eqn{( \gamma_0 t - \gamma_1 ) e^{-\gamma_2 t} + \gamma_1}
#'         over \[a, b), using a closed-form integral (with a small guard if \eqn{\gamma_2 \approx 0}).
#'   \item \code{"PiecewiseConstant"} — weighted average across overlapping
#'         piecewise intervals \[a_i, b_i) defined by \code{model_fixed_params$upper_cutoffs};
#'         one FOI per interval is taken from \code{par} (excluding any name \code{"rho", "k", "l", "w"}).
#'   \item \code{"Splines"} + \code{catalytic_model_type == "SimpleCatalytic"} —
#'         \eqn{\bar{\lambda}_{[a,b)} = \frac{1}{b-a} \int_a^b \frac{\pi'(t)}{1-\pi(t)} \, dt},
#'         computed by numerical integration applying \code{\link[stats]{predict}} to a
#'         fitted \code{smooth.spline}.
#'   \item \code{"Splines"} + \code{catalytic_model_type == "WaningImmunity"} —
#'         \eqn{\bar{\lambda}_{[a,b)} = \frac{1}{b-a} \int_a^b \frac{\pi'(t) + w\,\pi(t)}{1-\pi(t)} \, dt},
#'         where \code{w} is provided in \code{model_fixed_params$w}.
#' }
#'
#' If no specific form matches and \code{foi_t} is supplied, a generic numerical
#' mean is computed as \code{(1/(b-a)) * integrate(function(x) foi_t(x, par), a, b)}.
#' If \code{tau > 0}, the integral is taken over \[max(a, tau), b) and returns 0 when \code{b < tau}.
#'
#' @param catalytic_model_type Character. Model family (e.g., \code{"SimpleCatalytic"},
#'   \code{"WaningImmunity"}), used to disambiguate spline variants.
#' @param foi_functional_form Character. One of \code{"Constant"}, \code{"Linear"},
#'   \code{"Griffiths"}, \code{"Farringtons"}, \code{"PiecewiseConstant"}, \code{"Splines"}.
#' @param model_fixed_params Optional named list for required settings in some forms:
#'   \itemize{
#'     \item \code{Griffiths}: \code{list(tau = ...)}
#'     \item \code{PiecewiseConstant}: \code{list(upper_cutoffs = numeric())} — ascending upper bounds per interval
#'     \item \code{WaningImmunity + Splines}: \code{list(w = ...)}
#'   }
#' @param foi_t Optional function \code{foi_t(t, par)} used by the generic integration fallback
#'   when a closed-form or spline-specific method is not selected.
#' @param tau Numeric. Age threshold used for gating:
#'   \itemize{
#'     \item \code{Linear}, \code{Farringtons}: used only in the generic fallback branch.
#'     \item \code{Griffiths}: \code{model_fixed_params$tau} is used to gate FOI (returns 0 if \code{b <= tau}).
#'   }
#'
#' @return A function \code{group_foi(a, b, par)} for analytic forms, or
#'   \code{group_foi(a, b, spline_pi_t)} for spline forms, returning the mean FOI over \[a, b).
#'   If the combination is unrecognised and no \code{foi_t} is supplied, the helper will
#'   still return a function that numerically integrates \code{foi_t} if provided; otherwise it errors at call time.
#'
#' @details
#' For spline forms, \eqn{\pi(t)} and \eqn{\pi'(t)} are obtained via
#' \code{\link[stats]{predict}(spline_pi_t, t)} (with \code{deriv = 1} for the derivative), and
#' \eqn{\pi(t)} is capped below 1 to avoid division by zero. For \code{PiecewiseConstant},
#' the mean is a length-weighted average of interval FOIs by their overlap with \[a, b);
#' supply one FOI value per piecewise interval in \code{par} with the same order as \code{upper_cutoffs}.
#'
#' @note Helper for \code{\link{FoiFromCatalyticModel}}; not intended for independent use.
#'       Consider leaving this function unexported. If exported for advanced users,
#'       treat as internal API subject to change.
#'
#' @seealso \code{\link{set_foi_t}}, \code{\link{FoiFromCatalyticModel}}, \code{\link{FoiFromCatalyticModel_unparallelised}}
#'
#' @keywords internal
#' @importFrom stats integrate predict
#' @noRd
set_group_foi <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, foi_t = NULL, tau = 0) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (!is.na(foi_functional_form)) {
    if (foi_functional_form == "Constant") {
      group_foi <- function(a, b, par) {
        foi <- par[["foi"]]
        return(foi) # constant foi
      }
    }

    if (foi_functional_form == "Linear") {
      group_foi <- function(a, b, par) {
        m <- par[["m"]]
        c <- par[["c"]]

        return((m*a + c + m*b + c)/2)
      }
    }


    else if (foi_functional_form == "Griffiths") {
      tau <- model_fixed_params$tau
      group_foi <- function(a, b, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]

        if (b <= tau) return(0)

        # Adjust bounds if partially below tau
        a_adj <- max(a, tau)

        avg_foi<-gamma0*((b-a)/2 + a) + gamma0*gamma1
        return(avg_foi)
      }
    }

    else if (foi_functional_form == "Farringtons") {
      group_foi <- function(a, b, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        gamma2 <- par[["gamma2"]]
        if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
        integral <- (((gamma1 - b * gamma0) * gamma2 - gamma0) * exp(-b * gamma2) + exp(-a * gamma2) * ((b - a) * gamma1 * gamma2^2 * exp(a * gamma2) + (a * gamma0 - gamma1) * gamma2 + gamma0)) / gamma2^2
        avg_foi <- integral / (b - a)
        return(avg_foi)
      }
    }

    else if (foi_functional_form == "PiecewiseConstant") {
      upper_cutoffs <- model_fixed_params$upper_cutoffs
      lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
      group_foi <- function(a, b, par) {
        par <- par[ !(names(par) %in% c("rho", "k", "l", "w")) ]
        # Determine overlap between each piecewise interval and [a, b)
        overlaps <- pmin(b, upper_cutoffs) - pmax(a, lower_cutoffs)
        overlaps[overlaps < 0] <- 0  # No overlap = set to zero
        total_length <- b - a
        weighted_avg <- sum(par * overlaps) / total_length
        return(weighted_avg)
      }
    }

    else if (!is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Splines") {
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

    else if (!is.na(catalytic_model_type) && catalytic_model_type == "WaningImmunity" && foi_functional_form == "Splines") {
      w <- model_fixed_params$w
      group_foi <- function(a, b, spline_pi_t) {
        integrand <- function(x) {
          pi_x <- predict(spline_pi_t, x)$y
          dpi_x <- predict(spline_pi_t, x, deriv = 1)$y
          pi_x <- pmin(pi_x, 1 - 1e-8) # to prevent a divide by 0 error
          (dpi_x + pi_x*w) / (1 - pi_x)
        }
        result <- tryCatch({
          integrate(integrand, a, b)$value / (b - a) #  approximate the mean FOI over [a, b)  by integrating the FOI (π′ / (1−π)) over that interval and dividing by the interval width
        }, error = function(e) NA)
        return(result)
      }
    }

  }

  else if (tau > 0) {
    group_foi <- function(a, b, par) {
      if (b < tau) {
        return(0)
      }
      else if (a < tau) {
        1/(b-a) * integrate(function(x) {foi_t(x, par)}, tau, b)$value
      } else {
        1/(b-a) * integrate(function(x) {foi_t(x, par)}, a, b)$value
      }
    }
  }

  else {
    group_foi <- function(a, b, par) {
      1/(b-a) * integrate(function(x) {foi_t(x, par)}, a, b)$value
    }
  }

  return(group_foi)
}

