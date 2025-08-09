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
set_par_init <- function(type, model_fixed_params = NA) {
  if (type == "MuenchGeneral") {
    return(c(k=0.5, l=1, foi=0.1))
  }

  else if (type == "MuenchRestricted") {
    return(c(foi=0.1))
  }

  else if (type == "Griffiths") {
    return(c(gamma0=0.1, gamma1=-5))
  }

  else if (type == "Farringtons") {
    return(c(gamma0=0.1, gamma1=1, gamma2=0.1))
  }

  else if (type == "PiecewiseConstant") {
    num_pieces <- length(model_fixed_params$upper_cutoffs)
    par_vals <- rep(0.1, num_pieces)
    names(par_vals) <- paste0("foi", seq_len(num_pieces))
    return(par_vals)
  }

  else {
    return(NA)
  }
}
