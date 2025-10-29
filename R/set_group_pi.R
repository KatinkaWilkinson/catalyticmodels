#' Create a group prevalence function over \[a, b) (helper)
#'
#' Returns a function that computes the mean prevalence \eqn{\bar\pi_{[a,b)}} over
#' an age interval \[a, b) for a chosen catalytic-model form. This helper is used
#' internally by \code{\link{FoiFromCatalyticModel}} and related wrappers; it is
#' not intended for direct end-user calls.
#'
#' Supported combinations:
#' \itemize{
#'   \item \code{OriginalCatalytic + Constant} (with \code{tau == 0}) — closed form using
#'         \code{par[["k"]]}, \code{par[["l"]]}, \code{par[["foi"]]}.
#'   \item \code{RestrictedCatalytic + Constant} — uses fixed \code{k}, \code{l} from
#'         \code{model_fixed_params} and estimates \code{foi}; supports simple \code{tau} gating.
#'   \item \code{SimpleCatalytic + Constant} (with \code{tau == 0}) — closed form for mean prevalence.
#'   \item \code{SimpleCatalytic_NegativeCorrected + Linear} (with \code{tau == 0}) —
#'         integrates \code{1 - exp(-m t^2/2 - c t)}; uses \code{stats::pnorm} for \emph{erf}
#'         and \code{pracma::erfi} for \emph{erfi}. Returns 0 when the implied FOI is negative
#'         across \[a, b); splits at the root if the sign changes inside the interval.
#'   \item \code{SimpleCatalytic + Farringtons} (with \code{tau == 0}) — integrates a closed-form
#'         \eqn{\pi(t)} implied by the FOI \eqn{( \gamma_0 t - \gamma_1 ) e^{-\gamma_2 t} + \gamma_1};
#'         guarded when \eqn{\gamma_2 \approx 0}. On numeric failure returns \code{NA}.
#'   \item \code{SimpleCatalytic + PiecewiseConstant} — length-weighted average across overlapping
#'         piecewise intervals \[a_i, b_i) from \code{model_fixed_params$upper_cutoffs}. Parameter
#'         entries named \code{"rho"}, \code{"k"}, \code{"l"}, \code{"w"} are ignored when forming the
#'         per-piece FOIs.
#' }
#'
#' Fallback: if no specific form matches and \code{pi_t} is provided, the helper
#' returns a function that computes \code{(1/(b-a)) * integrate(function(x) pi_t(x, par), a, b)}.
#'
#' @param catalytic_model_type Character. Model family (e.g., \code{"OriginalCatalytic"},
#'   \code{"RestrictedCatalytic"}, \code{"SimpleCatalytic"}, \code{"SimpleCatalytic_NegativeCorrected"}).
#' @param foi_functional_form Character. One of \code{"Constant"}, \code{"Linear"},
#'   \code{"Farringtons"}, \code{"PiecewiseConstant"}.
#' @param model_fixed_params Optional named list of fixed settings required by some forms:
#'   \itemize{
#'     \item \code{RestrictedCatalytic}: \code{list(k = ..., l = ...)}
#'     \item \code{PiecewiseConstant}: \code{list(upper_cutoffs = numeric())} — ascending upper bounds per piece
#'   }
#' @param pi_t Optional function \code{pi_t(t, par)} used by the generic integration fallback
#'   when no closed-form/specific method is selected.
#' @param tau Numeric. Age threshold for gating. Behavior is model-specific in this helper:
#'   \itemize{
#'     \item In the \code{RestrictedCatalytic + Constant} branch, intervals entirely below \code{tau} return 0.
#'     \item In the \code{PiecewiseConstant} and generic fallback branches, the segment \[a, \min(b,\tau)) is
#'           treated as having \code{pi(t) = 1}, so the average over \[a, b) includes a contribution
#'           of \code{(\tau - a)} when \code{a < \tau < b}.
#'   }
#'
#' @return A function \code{group_pi(a, b, par)} that returns the mean prevalence over \[a, b).
#'   For piecewise models, \code{par} should contain one FOI value per piecewise interval (names optional),
#'   with entries named \code{"rho"}, \code{"k"}, \code{"l"}, \code{"w"} ignored in the computation.
#'   On numeric integration errors, \code{NA} may be returned.
#'
#' @details
#' \itemize{
#'   \item \strong{PiecewiseConstant}: the mean is a length-weighted average over overlaps
#'         with \[a, b). Overlap lengths are computed with \code{pmax/pmin}; negatives are clamped to 0.
#'   \item \strong{Farringtons}: mean prevalence is obtained by integrating a closed-form \eqn{\pi(t)};
#'         implemented with \code{\link[stats]{integrate}} and guarded for small \eqn{\gamma_2}.
#'   \item \strong{Linear (Negative-corrected)}: uses \emph{erf}/\emph{erfi} primitives; requires
#'         \pkg{pracma} for \code{erfi}. The result is 0 if the implied FOI is negative throughout \[a, b).
#' }
#'
#' @note Helper for \code{\link{FoiFromCatalyticModel}}; not intended for independent use.
#'       Consider leaving this function unexported. If exported for advanced users,
#'       treat it as internal API subject to change.
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}, \code{\link{set_group_foi}}, \code{\link{set_foi_t}}
#'
#' @keywords internal
#' @importFrom stats pnorm integrate
#' @noRd
set_group_pi <- function(catalytic_model_type, foi_functional_form, model_fixed_params, pi_t = NULL, tau=0) { # type is a string, model_fixed_params is a list
    if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "OriginalCatalytic" && foi_functional_form == "Constant" && tau == 0) {
      group_pi <- function(a, b, par) {
        k <- par[["k"]]
        l <- par[["l"]]
        foi <- par[["foi"]]
        if (foi < 1e-6) foi <- 1e-6
        return(k*(l-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a))))
      }
    }

    else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "RestrictedCatalytic" && foi_functional_form == "Constant" && tau == 0) {
      k <- model_fixed_params$k
      l <- model_fixed_params$l
      if (tau > 0) {
        group_pi <- function(a, b, par) {
          foi <- par[["foi"]]
          if (b < tau) {
            return(0)
          } else if (a < tau) {
            if (foi < 1e-6) foi <- 1e-6
            return(k*(l-(exp(-foi * tau) - exp(-foi * b)) / (foi * (b - a))))
          } else {
            if (foi < 1e-6) foi <- 1e-6
            return(k*(l-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a))))
          }
        }
      } else {
        group_pi <- function(a, b, par) {
          foi <- par[["foi"]]
          return(k*(l-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a))))
        }
      }
    }

  else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Constant" && tau == 0) {
    group_pi <- function(a, b, par) {
      foi <- par[["foi"]]
      return(1-(exp(-foi * a) - exp(-foi * b)) / (foi * (b - a)))
      # return(k * ((exp(-b * foi) * (b * l * foi * exp(b * foi) + 1)) / foi - (exp(-a * foi) * (a * l * foi * exp(a * foi) + 1)) / foi)) / (b - a)
    }
}

  else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic_NegativeCorrected" && foi_functional_form == "Linear" && tau == 0) {
    group_pi <- function(a, b, par) {
      m <- par[["m"]]
      c <- par[["c"]]

      # real erf via pnorm (no extra packages)
      erf <- function(x) {
        2 * pnorm(x * sqrt(2)) - 1
      }

      # real erfi via Dawson's integral: erfi(x) = 2/sqrt(pi) * exp(x^2) * dawson(x)
      # requires pracma for dawson(); avoids complex usage entirely
      erfi <- function(x) {
        pracma::erfi(x)
      }

      # --- main integral ---
      # Integrates 1 - exp(-m/2 * t^2 - c * t) from a to b
      integral_val <- function(a, b, c, m) {
        if (m > 0) {
          # m > 0: erf version
          s <- sqrt(m / 2)
          shift_b <- b + c / m
          shift_a <- a + c / m
          (b - a) -
            exp(c^2 / (2 * m)) * sqrt(pi / (2 * m)) *
            ( erf(s * shift_b) - erf(s * shift_a) )

        } else if (m < 0) {
          # m < 0: erfi version (all arguments real)
          s <- sqrt(-m / 2)
          shift_b <- b + c / m
          shift_a <- a + c / m
          (b - a) -
            exp(c^2 / (2 * m)) * sqrt(pi / (-2 * m)) *
            ( erfi(s * shift_b) - erfi(s * shift_a) )

        } else {
          # m == 0: integrand is 1 - exp(-c t)
          if (c != 0) {
            (b - a) + (exp(-c * b) - exp(-c * a)) / c
          } else {
            # m = 0, c = 0: integrand is 1 - exp(0) = 0
            0
          }
        }
      }

      if (m == 0) {
        if (c > 0) {
          return( integral_val(a, b, c, m) / (b-a) )
        } else {
          return(0)
        }
      }
      root <- -c/m
      if (a < root && b > root) { # then root is between a and b
        if (m > 0) { # positive slope therefore positive foi to the right of the root
          return( integral_val(root,b,c,m)/(b-a) )
        }
        else { # negative slope therefore positive foi to the left of the root
          return( integral_val(a,root,c,m)/(b-a) )
        }
      } else {
        if ((m*((a+b)/2) + c)>0) { # foi is positive in interval a,b
          return( integral_val(a,b,c,m)/(b-a) )
        }
        else (return(0)) # foi is negative in interval a,b
      }
    }
  }

    # else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Griffiths") {
    #   tau <- model_fixed_params$tau
    #   group_pi <- function(a, b, par) {
    #     gamma0 <- par[["gamma0"]]
    #     gamma1 <- par[["gamma1"]]
    #     m <- gamma0
    #     c <- gamma0*gamma1
    #
    #     # real erf via pnorm (no extra packages)
    #     erf <- function(x) {
    #       2 * pnorm(x * sqrt(2)) - 1
    #     }
    #
    #     # real erfi via Dawson's integral: erfi(x) = 2/sqrt(pi) * exp(x^2) * dawson(x)
    #     # requires pracma for dawson(); avoids complex usage entirely
    #     erfi <- function(x) {
    #       pracma::erfi(x)
    #     }
    #
    #     # --- main integral ---
    #     # Integrates 1 - exp(-m/2 * t^2 - c * t) from a to b
    #     integral_val <- function(a, b, c, m) {
    #       if (m > 0) {
    #         # m > 0: erf version
    #         s <- sqrt(m / 2)
    #         shift_b <- b + c / m
    #         shift_a <- a + c / m
    #         (b - a) -
    #           exp(c^2 / (2 * m)) * sqrt(pi / (2 * m)) *
    #           ( erf(s * shift_b) - erf(s * shift_a) )
    #
    #       } else if (m < 0) {
    #         # m < 0: erfi version (all arguments real)
    #         s <- sqrt(-m / 2)
    #         shift_b <- b + c / m
    #         shift_a <- a + c / m
    #         (b - a) -
    #           exp(c^2 / (2 * m)) * sqrt(pi / (-2 * m)) *
    #           ( erfi(s * shift_b) - erfi(s * shift_a) )
    #
    #       } else {
    #         # m == 0: integrand is 1 - exp(-c t)
    #         if (c != 0) {
    #           (b - a) + (exp(-c * b) - exp(-c * a)) / c
    #         } else {
    #           # m = 0, c = 0: integrand is 1 - exp(0) = 0
    #           0
    #         }
    #       }
    #     }
    #
    #     if (m == 0) {
    #       if (c > 0) {
    #         return( integral_val(a, b, c, m) / (b-a) )
    #       } else {
    #         return(0)
    #       }
    #     }
    #     root <- -c/m
    #     if (a < root && b > root) { # then root is between a and b
    #       if (m > 0) { # positive slope therefore positive foi to the right of the root
    #         return( integral_val(root,b,c,m)/(b-a) )
    #       }
    #       else { # negative slope therefore positive foi to the left of the root
    #         return( integral_val(a,root,c,m)/(b-a) )
    #       }
    #     } else {
    #       if ((m*((a+b)/2) + c)>0) { # foi is positive in interval a,b
    #         return( integral_val(a,b,c,m)/(b-a) )
    #       }
    #       else (return(0)) # foi is negative in interval a,b
    #     }
    #   }
    # }

    else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Farringtons" && tau == 0) {
      group_pi <- function(a, b, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        gamma2 <- par[["gamma2"]]
        if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
        pi_func <- function(t) {
          if (tau>0) {
            if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
            integral <- (((gamma1 - t * gamma0) * gamma2 - gamma0) * exp(-t * gamma2) + exp(-tau * gamma2) * ((t - tau) * gamma1 * gamma2^2 * exp(tau * gamma2) + (tau * gamma0 - gamma1) * gamma2 + gamma0)) / gamma2^2
            return(ifelse(t<tau, 1, 1-exp(-integral)))
          } else {
            if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
            return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
          }
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

    else if (!is.na(foi_functional_form) && !is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "PiecewiseConstant") {
      # upper_cutoffs <- model_fixed_params$upper_cutoffs
      # lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
      upper_cutoffs <- model_fixed_params$upper_cutoffs
      lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
      group_pi <- function(a, b, par) {
        # foi_pieces <- par[ names(par) != "rho" ]
        # # warning("Current foi_pieces: ", paste(round(foi_pieces, 4), collapse = ", "))
        # interval_lengths_a <- pmax(pmin(a, upper_cutoffs) - lower_cutoffs, 0)
        # cum_foi_a <- sum(foi_pieces * interval_lengths_a)
        # pi_a <- 1-exp(-cum_foi_a)
        #
        # interval_lengths_b <- pmax(pmin(b, upper_cutoffs) - lower_cutoffs, 0)
        # cum_foi_b <- sum(foi_pieces * interval_lengths_b)
        # pi_b <- 1-exp(-cum_foi_b)
        #

        integrand <- function(t, par) {
          foi_pieces <- unname(par[ !(names(par) %in% c("rho", "k", "l", "w")) ])
          sapply(t, function(tt) {
            interval_lengths <- pmax(pmin(tt, upper_cutoffs) - pmax(lower_cutoffs, tau), 0)
            cum_foi <- sum(foi_pieces * interval_lengths)
            1 - exp(-cum_foi)
          })
        }

        if (b <= tau) {
          return(1)
        } else if (a <= tau) {
          auc <- 1*(tau-a) + integrate(function(t) integrand(t, par), lower = tau, upper = b)$value

        } else {
          auc <- integrate(function(t) integrand(t, par), lower = a, upper = b)$value

        }
        return(auc/(b-a))
      }
    }

  else if (tau > 0) {
    group_pi <- function(a, b, par) {
      if (b <= tau) {
        return(1)
      } else if (a < tau) {
        1/(b-a) * (1*(tau-a) + integrate(function(x) {pi_t(x, par)}, tau, b)$value)
      } else {
        1/(b-a) * integrate(function(x) {pi_t(x, par)}, a, b)$value
      }
    }
  }


  else {
    group_pi <- function(a, b, par) {
      1/(b-a) * integrate(function(x) {pi_t(x, par)}, a, b)$value
    }
  }

  return(group_pi)
}
