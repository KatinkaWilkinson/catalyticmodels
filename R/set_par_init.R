#' Set initial parameter values for catalytic models (helper)
#'
#' Chooses default starting values for parameters based on the requested
#' \code{foi_functional_form} and \code{catalytic_model_type}. These are
#' pragmatic, non-zero seeds intended to help optimization routines converge.
#' This is a helper used internally by \code{FoiFromCatalyticModel()} and is
#' not typically called directly.
#'
#' @param catalytic_model_type Character or \code{NA}. High-level model family
#'   (e.g., \code{"OriginalCatalytic"}, \code{"RestrictedCatalytic"},
#'   \code{"SimpleCatalytic"}, \code{"WaningImmunity"}).
#' @param foi_functional_form Character or \code{NA}. FOI functional form used
#'   to pick sensible seeds. Supported here:
#'   \itemize{
#'     \item \code{"Constant"}: \code{c(foi = 0.1)}
#'     \item \code{"Linear"}: \code{c(m = -0.002, c = 0.12)}
#'     \item \code{"Griffiths"}: \code{c(gamma0 = -0.1, gamma1 = -1)}
#'     \item \code{"Farringtons"}: \code{c(gamma0 = 0.1, gamma1 = 1, gamma2 = 0.1)}
#'     \item \code{"PiecewiseConstant"}: length equals number of intervals, all \code{0.1}
#'           and named \code{foi1, foi2, ...}
#'   }
#' @param model_fixed_params Optional list of fixed values required by some
#'   configurations. For example:
#'   \itemize{
#'     \item \code{PiecewiseConstant}: \code{upper_cutoffs} (numeric vector) to
#'           determine how many \code{foi} pieces to seed.
#'     \item \code{WaningImmunity}: \code{w} (if known); if omitted or \code{NA},
#'           a default is added (see Details).
#'   }
#' @param rho Numeric scalar or \code{NA}. If \code{NA}, an initial value
#'   \code{rho = 0.8} is appended to the seed vector (and the caller should
#'   also expand parameter bounds accordingly).
#' @param pars Optional named numeric vector. A placeholder for externally
#'   supplied seeds. Note: for FOI-related parameters this helper constructs
#'   fresh defaults based on \code{foi_functional_form}; additional parameters
#'   specific to the catalytic family (e.g., \code{k}, \code{l}, \code{w},
#'   \code{rho}) are then appended as needed.
#'
#' @return A named numeric vector of initial parameter values. In addition to the
#'   FOI-related seeds listed above, the helper may append:
#'   \itemize{
#'     \item \code{c(k = 0.5, l = 0.5)} when \code{catalytic_model_type == "OriginalCatalytic"}.
#'     \item \code{c(w = 1/15)} when \code{catalytic_model_type == "WaningImmunity"} and
#'           \code{model_fixed_params$w} is missing or \code{NA}.
#'     \item \code{c(rho = 0.8)} when \code{rho} input is \code{NA}.
#'   }
#'
#' @details
#' \itemize{
#'   \item For \code{"PiecewiseConstant"}, the number of FOI seeds equals
#'         \code{length(model_fixed_params$upper_cutoffs)}; ensure this is provided.
#'   \item The chosen values are heuristics. You can override them upstream by
#'         providing your own \code{par_init} to \code{FoiFromCatalyticModel()}.
#'   \item This helper assembles a coherent starting vector that matches the
#'         expected parameterization of downstream likelihoods and link functions.
#' }
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}, \code{\link{set_foi_t}},
#'   \code{\link{set_pi_t}}, \code{\link{set_group_foi}}, \code{\link{set_group_pi}}
#' @noRd
set_par_init <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, rho, pars = NA) {
  if (!is.na(foi_functional_form)) {
    if (foi_functional_form == "Constant") {
      pars <- c(foi=0.1)
    }

    if (foi_functional_form == "Linear") {
      pars <- c(m=-0.002, c=0.12)
    }

    else if (foi_functional_form == "Griffiths") {
      pars <- c(gamma0=-0.1, gamma1=-1)
    }

    else if (foi_functional_form == "Farringtons") {
      pars <- c(gamma0=0.1, gamma1=1, gamma2=0.1)
    }

    else if (foi_functional_form == "PiecewiseConstant") {
      num_pieces <- length(model_fixed_params$upper_cutoffs)
      pars <- rep(0.1, num_pieces)
      names(pars) <- paste0("foi", seq_len(num_pieces))
    }
  }

  if (!is.na(catalytic_model_type)) {
    if (catalytic_model_type == "OriginalCatalytic") {
      pars <- c(pars, k=0.5, l=0.5)
    }
  }

  if (!is.na(catalytic_model_type)) {
    if (catalytic_model_type == "WaningImmunity") {
      if (is.na(model_fixed_params$w)) {
        pars <- c(pars, w=1/15)
      }
    }
  }

  if (is.na(rho)) {
    pars <- c(pars, rho=0.8)
  }

  return(pars)
}
