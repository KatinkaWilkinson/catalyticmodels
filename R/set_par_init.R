#' Set Initial Parameter Values for Catalytic Models
#'
#' Provides default starting values for model parameters based on the specified catalytic model type.
#' These initial values are useful for optimization routines, such as maximum likelihood estimation,
#' to ensure reasonable starting points for parameter searches.
#'
#' @param type Character string specifying the model type. Supported values include:
#'   \code{"MuenchGeneral"}, \code{"MuenchRestricted"}, \code{"Griffiths"},
#'   \code{"Farringtons"}, and \code{"PiecewiseConstant"}.
#' @param model_fixed_params Optional named list of fixed parameters required for certain types:
#'   for example, \code{"PiecewiseConstant"} requires \code{upper_cutoffs} to determine the number of intervals.
#'
#' @return A named numeric vector of initial parameter values:
#' \itemize{
#'   \item \code{"MuenchGeneral"}: \code{c(k=0.5, l=1, foi=0.1)}
#'   \item \code{"MuenchRestricted"}: \code{c(foi=0.1)}
#'   \item \code{"Griffiths"}: \code{c(gamma0=0.1, gamma1=-5)}
#'   \item \code{"Farringtons"}: \code{c(gamma0=0.1, gamma1=1, gamma2=0.1)}
#'   \item \code{"PiecewiseConstant"}: A vector of 0.1 for each interval, named \eqn{foi1, foi2, \dots}
#' }
#' If \code{type} is not recognized, the function returns \code{NA}.
#'
#' @details
#' For \code{"PiecewiseConstant"}, the number of parameters is determined by the length of
#' \code{model_fixed_params\$upper_cutoffs}. Ensure this list element exists to avoid errors.
#' The default values are chosen to be moderate and non-zero, aiding in stable optimisation.
#' @export
set_par_init <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, rho, pars = NA) {
  if (!is.na(foi_functional_form)) {
    if (foi_functional_form == "Constant") {
      pars <- c(foi=0.1)
    }

    else if (foi_functional_form == "Griffiths") {
      pars <- c(gamma0=1, gamma1=-1)
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

  if (is.na(rho)) {
    pars <- c(pars, rho=0.8) # a reasonable assay sensitivity, based on the results that popped up when I googled what common assay sensitivities were
  }

  return(pars)
}
