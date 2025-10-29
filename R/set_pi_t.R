#' Create a prevalence function \eqn{\pi(t)} for a chosen catalytic form (helper)
#'
#' Returns a function \code{pi_t(t, par)} that computes modelled prevalence at
#' age \code{t} for the specified \code{catalytic_model_type} and
#' \code{foi_functional_form}. This is a helper used internally by
#' \code{\link{FoiFromCatalyticModel}}; it is not intended for direct end-user calls.
#'
#' Supported combinations (as implemented here):
#' \itemize{
#'   \item \code{OriginalCatalytic + Constant}: \eqn{\pi(t)=k\{\,l-\exp(-\lambda t)\,\}}, uses \code{par[["k"]]}, \code{par[["l"]]}, \code{par[["foi"]]}.
#'   \item \code{RestrictedCatalytic + Constant}: same form with fixed \code{k,l} from \code{model_fixed_params}.
#'   \item \code{SimpleCatalytic + Griffiths}: linear FOI after \eqn{\tau}, \eqn{\lambda(t)=\gamma_0(t+\gamma_1)}; integrates piecewise.
#'   \item \code{SimpleCatalytic + Farringtons}: damped linear FOI
#'         \eqn{(\gamma_0 t-\gamma_1)\exp(-\gamma_2 t)+\gamma_1} with a closed-form \eqn{\pi(t)}.
#'   \item \code{SimpleCatalytic + PiecewiseConstant}: piecewise \eqn{\lambda} over
#'         \code{model_fixed_params$upper_cutoffs}; expects one FOI per piece in \code{par}
#'         (entries named \code{"rho"} are ignored).
#'   \item \code{SimpleCatalytic_NegativeCorrected + Linear}: FOI \eqn{\lambda(t)=mt+c} with
#'         sign handling; splits at zero crossings (roots) found via
#'         \code{rootSolve::uniroot.all} and accumulates only positive-FOI segments
#'         using \code{group_foi}.
#'   \item \code{OriginalCatalytic}, \code{SimpleCatalytic}, \code{RestrictedCatalytic}
#'         (generic fallback): \eqn{\pi(t)=1-\exp\{-\int_0^t \lambda(x)\,dx\}} or its
#'         \code{k,l} analogue, evaluated with \code{stats::integrate}.
#'   \item \code{WaningImmunity} (and \code{WaningImmunity_NegativeCorrected}):
#'         \eqn{\pi(t)=1-\big[\exp(-\int_0^t(\lambda+w)) + w\int_0^t \exp(-\int_u^t(\lambda+w))\,du\big]},
#'         where \code{w} is provided in \code{model_fixed_params$w} or as \code{par[["w"]]}.
#' }
#'
#' \strong{Note on \eqn{\tau} gating:} In all branches here, if \code{t < tau} the function
#' returns \code{1} (maternal antibodies), i.e. \eqn{\pi(t)=1} for \eqn{t<\tau}.
#'
#' @param catalytic_model_type Character. Model family (e.g., \code{"OriginalCatalytic"},
#'   \code{"RestrictedCatalytic"}, \code{"SimpleCatalytic"},
#'   \code{"SimpleCatalytic_NegativeCorrected"}, \code{"WaningImmunity"},
#'   \code{"WaningImmunity_NegativeCorrected"}).
#' @param foi_functional_form Character. One of \code{"Constant"}, \code{"Linear"},
#'   \code{"Griffiths"}, \code{"Farringtons"}, \code{"PiecewiseConstant"}.
#'   (Spline-based prevalence is handled elsewhere; this helper does not create a spline \code{pi_t}.)
#' @param model_fixed_params Optional named list of fixed values required by some forms:
#'   \itemize{
#'     \item \code{RestrictedCatalytic}: \code{list(k = ..., l = ...)}
#'     \item \code{Griffiths}: \code{list(tau = ...)}
#'     \item \code{PiecewiseConstant}: \code{list(upper_cutoffs = numeric())}
#'     \item \code{WaningImmunity*}: \code{list(w = ...)} if \code{w} is known
#'   }
#' @param foi_t Optional FOI function \code{foi_t(t, par)} used by the fallback
#'   integral forms and by the negative-corrected linear variant.
#' @param tau Numeric. Maternal-antibody threshold; this helper returns \code{1} for
#'   \code{t < tau} in all branches.
#'
#' @return A function \code{pi_t(t, par)} returning prevalence for each \code{t}.
#'   For piecewise forms, \code{par} should contain one FOI value per piece in the same
#'   order as \code{upper_cutoffs} (names optional). Entries named \code{"rho"} are ignored.
#'   For waning-immunity forms, \code{w} must be supplied in \code{model_fixed_params}
#'   or as \code{par[["w"]]}.
#'
#' @details
#' \itemize{
#'   \item \strong{Integral forms}: use \code{stats::integrate}; numerical tolerances are left at defaults.
#'   \item \strong{Negative-corrected linear}: roots are found on \code{c(0, t]} via
#'         \code{rootSolve::uniroot.all}; only segments where \eqn{\lambda(t)>0} contribute.
#'   \item \strong{PiecewiseConstant}: prevalence is \eqn{1-\exp(-\sum_i \lambda_i L_i(t))}, where
#'         \eqn{L_i(t)} is the overlap length of \code{t} with piece \code{i} after applying \code{tau}.
#' }
#'
#' @note Helper for \code{\link{FoiFromCatalyticModel}}; consider leaving it unexported.
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}, \code{\link{set_group_pi}},
#'   \code{\link{set_foi_t}}
#'
#' @keywords internal
#' @importFrom stats integrate
#' @importFrom rootSolve uniroot.all
#' @noRd
set_pi_t <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, foi_t=NULL, tau=0) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (!is.na(foi_functional_form) && catalytic_model_type == "OriginalCatalytic" && foi_functional_form == "Constant") {
    pi_t <- function(t, par) {
      k <- par[["k"]]
      l <- par[["l"]]
      foi <- par[["foi"]]
      return(ifelse(t<tau,1,k * (l - exp(-foi * t))))
    }
  }

  else if (!is.na(foi_functional_form) && catalytic_model_type == "RestrictedCatalytic" && foi_functional_form == "Constant") {
    k <- model_fixed_params$k
    l <- model_fixed_params$l
    pi_t <- function(t, par) {
      foi <- par[["foi"]]
      return(ifelse(t<tau,1,k * (l - exp(-foi * t))))
    }
  }

  else if (!is.na(foi_functional_form) && catalytic_model_type == "SimpleCatalytic_NegativeCorrected" && foi_functional_form == "Linear") { # linear with maternal antibodies is the actual form...
    pi_t <- function(t, par) {
      pi_t_one <- function(tt) {
        if (tt < tau) { return(1) }
        m <- par[["m"]]
        c <- par[["c"]]

        mean_foi <- function(a, b) {
          if (b < tau) {
            return(0)
          } else if (a < tau) {
            (m*tau + c + m*b + c)/2 * (b-tau)/(b-a)
          } else {
            (m*a + c + m*b + c)/2
          }
        }

        if (m == 0) {
          if (c > 0) {
            return(1 - exp(-c * tt))
          } else {
            return(0)
          }
        }
        root <- -c / m
        if (root < tt) {
          if (foi_t(tt, par) > 0) {
            return(1 - exp(-(tt - root) * mean_foi(root, tt)))
          } else {
            return(1 - exp(-(root) * mean_foi(0, root)))
          }
        } else {
          if (foi_t(tt, par) > 0) {
            return(1 - exp(-(tt) * mean_foi(0, tt)))
          } else {
            return(0)
          }
        }
      }

      if (length(t) == 1L) {
        pi_t_one(t)
      } else {
        vapply(t, pi_t_one, numeric(1))
      }
    }
  }

  else if (!is.na(foi_functional_form) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Griffiths") { # linear with maternal antibodies is the actual form...
    tau <- model_fixed_params$tau
    pi_t <- function(t, par) {
      pi_t_one <- function(tt) {
        if (tt < tau) { return(1) }
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        m <- gamma0
        c <- gamma0*gamma1


        mean_foi <- function(a, b) {
          if (b < tau) {
            return(0)
          } else if (a < tau) {
            (m*tau + c + m*b + c)/2 * (b-tau)/(b-a)
          } else {
            (m*a + c + m*b + c)/2
          }
        }

        if (m == 0) {
          if (c > 0) {
            return(1 - exp(-c * tt))
          } else {
            return(0)
          }
        }
        root <- -c / m
        if (root < tt) {
          if (foi_t(tt, par) > 0) {
            return(1 - exp(-(tt - root) * mean_foi(root, tt)))
          } else {
            return(1 - exp(-(root) * mean_foi(0, root)))
          }
        } else {
          if (foi_t(tt, par) > 0) {
            return(1 - exp(-(tt) * mean_foi(0, tt)))
          } else {
            return(0)
          }
        }
      }

      if (length(t) == 1L) {
        pi_t_one(t)
      } else {
        vapply(t, pi_t_one, numeric(1))
      }
    }

  }

  else if (!is.na(foi_functional_form) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Farringtons") { # linear dampered is the actual functional form
    pi_t <- function(t, par) {
      gamma0 <- par[["gamma0"]]
      gamma1 <- par[["gamma1"]]
      gamma2 <- par[["gamma2"]]
      if (tau>0) {
        if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
        integral <- (((gamma1 - t * gamma0) * gamma2 - gamma0) * exp(-t * gamma2) + exp(-tau * gamma2) * ((t - tau) * gamma1 * gamma2^2 * exp(tau * gamma2) + (tau * gamma0 - gamma1) * gamma2 + gamma0)) / gamma2^2
        return(ifelse(t<tau, 1, 1-exp(-integral)))
      } else {
        return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
      }
    }
  }

  else if (!is.na(foi_functional_form) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])

    pi_t <- function(t, par) {
      foi_pieces <- par[ names(par) != "rho" ]
      interval_lengths <- pmax(pmin(t, upper_cutoffs) - max(tau, lower_cutoffs), 0)
      cum_foi <- sum(foi_pieces * interval_lengths)
      return(ifelse(t<tau, 1, 1 - exp(-cum_foi)))
    }
  }

  else if (catalytic_model_type == "OriginalCatalytic") { # potentiallyyyy incorporate k and l... do they also apply to the other models, though?
    # pi_t: compute π(t) = 1 - exp( -∫_0^t λ(x) dx )
    pi_t <- function(t, par) {
      k <- par[["k"]]
      l <- par[["l"]]
      pi_one <- function(tt) {
        if (tt < tau) { return(1) }
        I0 <- stats::integrate(foi_t, lower = tau, upper = tt, par=par)$value
        val <- k*(l - exp(-I0))
        # max(0, min(1, val))  # clamp tiny numerical drift
      }
      vapply(t, pi_one, numeric(1))
    }
  }

  else if (catalytic_model_type == "SimpleCatalytic") {
    # pi_t: compute π(t) = 1 - exp( -∫_0^t λ(x) dx )
    pi_t <- function(t, par) {
      pi_one <- function(tt) {
        if (tt < tau) { return(1) }
        I0 <- stats::integrate(foi_t, lower = tau, upper = tt, par=par)$value
        val <- 1 - exp(-I0)
        # max(0, min(1, val))  # clamp tiny numerical drift
      }
      vapply(t, pi_one, numeric(1))
    }
  }

  # I must decide how to handle k and l! Will I have an original catalytic and a restricted one??
  else if (catalytic_model_type == "RestrictedCatalytic") {
    k <- model_fixed_params$k
    l <- model_fixed_params$l
    # pi_t: compute π(t) = 1 - exp( -∫_0^t λ(x) dx )
    pi_t <- function(t, par) {
      pi_one <- function(tt) {
        if (tt < tau) { return(1) }
        I0 <- stats::integrate(foi_t, lower = tau, upper = tt, par=par)$value
        val <- k*(l - exp(-I0))
        # max(0, min(1, val))  # clamp tiny numerical drift
      }
      vapply(t, pi_one, numeric(1))
    }
  }

  else if (catalytic_model_type == "SimpleCatalytic_NegativeCorrected") {
    pi_t <- function(t, par) {
      pi_one <- function(tt, par) {
        if (tt < tau) { return(1) }
        roots <- sort(rootSolve::uniroot.all(function(x) foi_t(x, par), c(0.00001, tt)))

        # Handle case with no roots
        if (length(roots) == 0) {
          if (foi_t(tt/2, par) > 0) {
            return((tt) * group_foi(0, tt, par))
          } else {
            return(0)
          }
        }

        cumulative_foi <- 0

        # First segment [0, root[1]]
        if (foi_t(roots[1] / 2, par) > 0) {
          cumulative_foi <- cumulative_foi + (roots[1])*group_foi(0, roots[1], par)
        }

        # Segments between roots
        for (i in seq_len(length(roots) - 1)) {
          mid <- (roots[i] + roots[i+1]) / 2
          if (foi_t(mid, par) > 0) {
            cumulative_foi <- cumulative_foi + (roots[i+1]-roots[i])*group_foi(roots[i], roots[i+1], par)
          }
        }

        # Final segment [last root, tt]
        if (foi_t((roots[length(roots)] + tt) / 2, par) > 0) {
          cumulative_foi <- cumulative_foi + (tt-roots[length(roots)])*group_foi(roots[length(roots)], tt, par)
        }

        val <- 1 - exp(-cumulative_foi)
        return(val)
      }
      vapply(t, pi_one, numeric(1))
      }
  }


  else if (catalytic_model_type == "WaningImmunity" && foi_functional_form == "PiecewiseConstant") {
    vec1 <- function(f) {
      function(x, ...) {
        y <- f(x, ...)
        if (length(y) == length(x)) return(y)
        if (length(y) == 1L)       return(rep_len(y, length(x)))
        vapply(x, function(xi) f(xi, ...), numeric(1))
      }
    }
    par <- par[ !(names(par) %in% c("rho", "k", "l", "w")) ]
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    w <- model_fixed_params$w  # seroreversion rate
    if (!is.na(w)) { # w is known!
      # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
      pi_t <- function(t, par, rel.tol = .Machine$double.eps^0.5) {
        # λ(x)+w
        integrand <- vec1(function(x) foi_t(x, par) + w)

        integrand <- function(x) {
          # Determine overlap between each piecewise interval and [a, b)
          overlaps <- pmin(b, upper_cutoffs) - pmax(a, lower_cutoffs)
          overlaps[overlaps < 0] <- 0  # No overlap = set to zero
          total_length <- b - a
          cumulative <- sum(par * overlaps) + total_length*w
        }

        pi_one <- function(tt) {
          if (tt < tau) { return(1) }

          if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

          # I0 = ∫_0^t (λ+w) dx
          I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

          # inner_exp(u)
          inner_exp <- function(u) {
            vapply(u, function(ui) {
              Iu <- stats::integrate(integrand, lower = ui, upper = tt, rel.tol = rel.tol)$value
              exp(-Iu)
            }, numeric(1))
          }

          # term2 = w * ∫_0^t exp(-∫_u^t (λ+w)) du
          term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt, rel.tol = rel.tol)$value

          val <- 1 - (exp(-I0) + term2)
          # optional clamp to [0,1]:
          # val <- max(0, min(1, val))
          return(val)
        }

        vapply(t, pi_one, numeric(1))
      }
    } else { # w = is unknown
      # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
      pi_t <- function(t, par, rel.tol = .Machine$double.eps^0.5) {
        w <- par[["w"]]
        # λ(x)+w
        integrand <- vec1(function(x) foi_t(x, par) + w)

        pi_one <- function(tt) {
          if (tt < tau) { return(1) }

          if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

          # I0 = ∫_0^t (λ+w) dx
          I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

          # inner_exp(u):
          inner_exp <- function(u) {
            vapply(u, function(ui) {
              Iu <- stats::integrate(integrand, lower = ui, upper = tt, rel.tol = rel.tol)$value
              exp(-Iu)
            }, numeric(1))
          }

          # term2 = w * ∫_0^t exp(-∫_u^t (λ+w)) du
          term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt, rel.tol = rel.tol)$value

          val <- 1 - (exp(-I0) + term2)
          # optional clamp to [0,1]:
          # val <- max(0, min(1, val))
          return(val)
        }

        vapply(t, pi_one, numeric(1))
      }
    }



  }

  else if (catalytic_model_type == "WaningImmunity") {
    vec1 <- function(f) {
      function(x, ...) {
        y <- f(x, ...)
        if (length(y) == length(x)) return(y)
        if (length(y) == 1L)       return(rep_len(y, length(x)))
        vapply(x, function(xi) f(xi, ...), numeric(1))
      }
    }

    w <- model_fixed_params$w  # seroreversion rate
    if (!is.na(w)) { # w is known!
      # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
      pi_t <- function(t, par, rel.tol = .Machine$double.eps^0.5) {
        # λ(x)+w
        integrand <- vec1(function(x) foi_t(x, par) + w)

        pi_one <- function(tt) {
          if (tt < tau) { return(1) }

          if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

          # I0 = ∫_0^t (λ+w) dx
          I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

          # inner_exp(u)
          inner_exp <- function(u) {
            vapply(u, function(ui) {
              Iu <- stats::integrate(integrand, lower = ui, upper = tt, rel.tol = rel.tol)$value
              exp(-Iu)
            }, numeric(1))
          }

          # term2 = w * ∫_0^t exp(-∫_u^t (λ+w)) du
          term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt, rel.tol = rel.tol)$value

          val <- 1 - (exp(-I0) + term2)
          # optional clamp to [0,1]:
          # val <- max(0, min(1, val))
          return(val)
        }

        vapply(t, pi_one, numeric(1))
      }
    } else { # w = is unknown
      # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
      pi_t <- function(t, par, rel.tol = .Machine$double.eps^0.5) {
        w <- par[["w"]]
        # λ(x)+w
        integrand <- vec1(function(x) foi_t(x, par) + w)

        pi_one <- function(tt) {
          if (tt < tau) { return(1) }

          if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

          # I0 = ∫_0^t (λ+w) dx
          I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

          # inner_exp(u):
          inner_exp <- function(u) {
            vapply(u, function(ui) {
              Iu <- stats::integrate(integrand, lower = ui, upper = tt, rel.tol = rel.tol)$value
              exp(-Iu)
            }, numeric(1))
          }

          # term2 = w * ∫_0^t exp(-∫_u^t (λ+w)) du
          term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt, rel.tol = rel.tol)$value

          val <- 1 - (exp(-I0) + term2)
          # optional clamp to [0,1]:
          # val <- max(0, min(1, val))
          return(val)
        }

        vapply(t, pi_one, numeric(1))
      }
    }



  }

  else if (catalytic_model_type == "WaningImmunity_NegativeCorrected") {
    w <- model_fixed_params$w  # seroreversion rate
    vec1 <- function(f) {
      function(x, ...) {
        y <- f(x, ...)
        if (length(y) == length(x)) return(y)
        if (length(y) == 1L)       return(rep_len(y, length(x)))
        vapply(x, function(xi) f(xi, ...), numeric(1))
      }
    }

    # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
    pi_t <- function(t, par, rel.tol = 1e-4) {
      # λ(x)+w
      integrand <- vec1(function(x) max(0,foi_t(x, par)) + w)

      pi_one <- function(tt) {
        if (tt < tau) { return(1) }

        if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

        # I0 = ∫_0^t (λ+w) dx
        I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

        # inner_exp(u)
        inner_exp <- function(u) {
          vapply(u, function(ui) {
            Iu <- stats::integrate(integrand, lower = ui, upper = tt, rel.tol = rel.tol)$value
            exp(-Iu)
          }, numeric(1))
        }

        # term2 = w * ∫_0^t exp(-∫_u^t (λ+w)) du
        term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt, rel.tol = rel.tol)$value

        val <- 1 - (exp(-I0) + term2)
        # optional clamp to [0,1]:
        # val <- max(0, min(1, val))
        return(val)
      }

      vapply(t, pi_one, numeric(1))
    }
  }


  else {
    return(NA)
  }

  return(pi_t)
}
