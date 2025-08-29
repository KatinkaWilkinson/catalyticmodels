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
set_pi_t <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, foi_t=NULL) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (!is.na(foi_functional_form)) {
    if (catalytic_model_type == "OriginalCatalytic" && foi_functional_form == "Constant") {
      pi_t <- function(t, par) {
        k <- par[["k"]]
        l <- par[["l"]]
        foi <- par[["foi"]]
        return(k * (l - exp(-foi * t)))
      }
    }

    else if (catalytic_model_type == "RestrictedCatalytic" && foi_functional_form == "Constant") {
      k <- model_fixed_params$k
      l <- model_fixed_params$l
      pi_t <- function(t, par) {
        foi <- par[["foi"]]
        return(k * (l - exp(-foi * t)))
      }
    }

    else if (catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Griffiths") { # linear with maternal antibodies is the actual form...
      tau <- model_fixed_params$tau
      pi_t <- function(t, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        return(ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau)))))
      }
    }

    else if (catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Farringtons") { # linear dampered is the actual functional form
      pi_t <- function(t, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        gamma2 <- par[["gamma2"]]
        return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
      }
    }

    else if (catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "PiecewiseConstant") {
      upper_cutoffs <- model_fixed_params$upper_cutoffs
      lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])

      pi_t <- function(t, par) {
        foi_pieces <- par[ names(par) != "rho" ]
        interval_lengths <- pmax(pmin(t, upper_cutoffs) - lower_cutoffs, 0)
        cum_foi <- sum(foi_pieces * interval_lengths)
        return(1 - exp(-cum_foi))
      }
    }
    return(pi_t)
  }


  else if (catalytic_model_type == "OriginalCatalytic") { # potentiallyyyy incorporate k and l... do they also apply to the other models, though?
    # pi_t: compute π(t) = 1 - exp( -∫_0^t λ(x) dx )
    pi_t <- function(t, par) {
      k <- par[["k"]]
      l <- par[["l"]]
      pi_one <- function(tt) {
        I0 <- stats::integrate(foi_t, lower = 0, upper = tt, par=par)$value
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
        I0 <- stats::integrate(foi_t, lower = 0, upper = tt, par=par)$value
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
        I0 <- stats::integrate(foi_t, lower = 0, upper = tt, par=par)$value
        val <- k*(l - exp(-I0))
        # max(0, min(1, val))  # clamp tiny numerical drift
      }
      vapply(t, pi_one, numeric(1))
    }
  }

  # else if (catalytic_model_type == "WaningImmunity") {
  #     w <- model_fixed_params$w
  #     # pi_t: compute π(t) = 1 - ( exp(-∫_0^t [λ(x)+w] dx) + w ∫_0^t exp(-∫_u^t [λ(x)+w] dx) du )
  #     pi_t <- function(t, par) {
  #       pi_one <- function(tt) {
  #         # integrand for ∫(λ(x)+w)dx
  #         integrand <- function(x) foi_t(x, par) + w
  #
  #         # I0 = ∫_0^t (λ(x)+w) dx
  #         I0 <- stats::integrate(integrand, lower = 0, upper = tt)$value
  #
  #         # inner: for each u, compute exp( - ∫_u^t (λ(x)+w) dx )
  #         inner_exp <- function(u) {
  #           Iu <- stats::integrate(integrand, lower = u, upper = tt)$value
  #           exp(-Iu)
  #         }
  #
  #         # term2 = w * ∫_0^t exp( - ∫_u^t (λ+ w) ) du
  #         term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt)$value
  #
  #         # π(t) = 1 - ( e^{-I0} + term2 )
  #         val <- 1 - (exp(-I0) + term2)
  #
  #         # clamp tiny numerical drift outside [0,1]
  #         #max(0, min(1, val))
  #       }
  #
  #     vapply(t, pi_one, numeric(1))
  #   }
  # }

  else if (catalytic_model_type == "WaningImmunity") {
    w <- model_fixed_params$w  # seroreversion rate

    # helper: force any 1D function to return a vector matching x
    vec1 <- function(f) {
      function(x, ...) {
        y <- f(x, ...)
        if (length(y) == length(x)) return(y)
        if (length(y) == 1L)       return(rep_len(y, length(x)))
        vapply(x, function(xi) f(xi, ...), numeric(1))
      }
    }

    # π(t) = 1 - ( exp(-∫0^t (λ+w)) + w ∫0^t exp(-∫u^t (λ+w)) du )
    pi_t <- function(t, par, rel.tol = .Machine$double.eps^0.5) {
      # vectorised λ(x)+w
      integrand <- vec1(function(x) foi_t(x, par) + w)

      pi_one <- function(tt) {
        if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")

        # I0 = ∫_0^t (λ+w) dx
        I0 <- stats::integrate(integrand, lower = 0, upper = tt, rel.tol = rel.tol)$value

        # inner_exp(u): return vector same length as u
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


  else if (catalytic_model_type == "Vaccine") {
    # ------------------------------------------------------------
    # Model 1 (continuous-hazard vaccination)
    # pi_t: compute π(t) under dπ/dt = (1-π)[λ(t) + ε ν(t)] - ω π
    # Closed form:
    #   Let Λ(t) = λ(t) + ε ν(t).
    #   π(t) = 1 - (1 - π0) * exp(-∫_0^t [Λ(x) + ω] dx)
    #                - ω * ∫_0^t exp(-∫_u^t [Λ(x) + ω] dx) du
    # Requirements:
    #   - foi_t: function of time t -> λ(t)
    #   - v_t:  function of time t -> ν(t)    (vaccination rate)
    #   - t:     numeric scalar or vector of times
    #   - eps:   vaccination seroconversion efficacy ε (default 1)
    #   - w:     waning rate ω >= 0             (default 0)
    #   - pi0:   initial seropositive fraction π(0) (default 0)
    # ------------------------------------------------------------
    v_t <- model_fixed_params$v_t
    eps <- model_fixed_params$eps
    w <- model_fixed_params$w
    pi0 <- model_fixed_params$pi0

    pi_t <- function(t, par) {

      if (!is.function(foi_t)) stop("foi_t must be a function of one numeric argument.")
      if (!is.function(v_t))  stop("v_t must be a function of one numeric argument.")
      stopifnot(is.numeric(eps), length(eps) == 1,
                is.numeric(w),   length(w)   == 1,
                is.numeric(pi0), length(pi0) == 1)

      # Λ(t) = λ(t) + ε ν(t)
      bigLambda <- function(x) foi_t(x, par) + eps * v_t(x)

      pi_one <- function(tt) {
        if (!is.finite(tt) || tt < 0) stop("t must be nonnegative and finite.")
        integrand <- function(x) bigLambda(x) + w

        # ∫_0^t [Λ(x) + ω] dx
        I0 <- stats::integrate(integrand, lower = 0, upper = tt)$value

        # inner: exp( - ∫_u^t [Λ(x) + ω] dx )
        inner_exp <- function(u) {
          Iu <- stats::integrate(integrand, lower = u, upper = tt)$value
          exp(-Iu)
        }

        term2 <- w * stats::integrate(inner_exp, lower = 0, upper = tt)$value
        val   <- 1 - (1 - pi0) * exp(-I0) - term2

        # clamp tiny numerical drift
        max(0, min(1, val))
      }

      vapply(t, pi_one, numeric(1))
    }

  }

  else {
    return(NA)
  }

  return(pi_t)
}
