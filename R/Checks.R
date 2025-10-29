#' Validate Inputs for Catalytic-Model Fitting (internal)
#'
#' Performs comprehensive argument validation for functions that fit catalytic
#' models (e.g., \code{FoiFromCatalyticModel}). On success it returns
#' \code{TRUE} invisibly; otherwise it throws an error explaining all detected
#' issues. Some non-fatal situations produce warnings.
#'
#' @param t Either a numeric vector of exact ages, or a numeric two-column
#'   matrix of age intervals \[a, b\].
#' @param y Numeric vector of seropositive counts, same length as the number of
#'   age rows implied by \code{t}.
#' @param n Numeric vector of total tested counts, same length as \code{y}.
#' @param pi_t Optional function \code{pi_t(t, par)} returning prevalence at age
#'   \code{t} for parameter vector \code{par}.
#' @param foi_t Optional function \code{foi_t(t, par)} returning instantaneous
#'   force of infection at age \code{t}.
#' @param group_pi Optional function \code{group_pi(a, b, par)} returning mean
#'   prevalence over an interval \[a, b\].
#' @param group_foi Optional function \code{group_foi(a, b, par)} returning mean
#'   force of infection over an interval \[a, b\].
#' @param par_init Optional named numeric vector of initial parameter values, required
#'   when any of \code{pi_t}, \code{group_pi}, \code{foi_t}, or \code{group_foi} is supplied.
#' @param rho Either \code{NA} (to estimate) or a numeric scalar in \[0, 1\] for
#'   test adjustment.
#' @param catalytic_model_type Optional character string selecting a preset model
#'   family. Allowed values: \code{"OriginalCatalytic"}, \code{"RestrictedCatalytic"},
#'   \code{"SimpleCatalytic"}, \code{"WaningImmunity"},
#'   \code{"SimpleCatalytic_NegativeCorrected"}, \code{"WaningImmunity_NegativeCorrected"}.
#' @param foi_functional_form Optional character string selecting a preset FOI
#'   form. Allowed values: \code{"Constant"}, \code{"Linear"}, \code{"Griffiths"},
#'   \code{"Farringtons"}, \code{"PiecewiseConstant"}, \code{"Splines"}.
#' @param model_fixed_params Optional list of fixed hyperparameters required by
#'   some presets, e.g.,
#'   \itemize{
#'     \item \code{Griffiths}: \code{list(tau = numeric_scalar)}
#'     \item \code{PiecewiseConstant}: \code{list(upper_cutoffs = numeric_vector)}
#'     \item \code{RestrictedCatalytic}: \code{list(k = scalar, l = scalar)}
#'     \item \code{WaningImmunity}: \code{list(w = scalar\ or\ NA)}
#'   }
#' @param lower,upper Lower and upper bounds used by optimisation. When using
#'   one-parameter Brent optimisation, both must be finite with \code{lower < upper}.
#'
#' @details
#' The following checks are performed:
#' \enumerate{
#'   \item \strong{Shapes}: \code{t} must be a numeric vector or a numeric
#'         two-column matrix. The lengths of \code{y} and \code{n} must match
#'         the number of rows implied by \code{t}.
#'   \item \strong{Model choice consistency}:
#'         \itemize{
#'           \item For vector \code{t}, supply exactly one of \code{pi_t} or
#'                 \code{catalytic_model_type} (not both).
#'           \item Supply exactly one of \code{foi_t} or
#'                 \code{foi_functional_form} (not both).
#'           \item Provide at least one of \code{pi_t}, \code{group_pi}, or
#'                 \code{catalytic_model_type}. If \code{catalytic_model_type}
#'                 is provided, do not also provide \code{pi_t} or
#'                 \code{group_pi}.
#'         }
#'   \item \strong{Function signatures}: If provided, \code{pi_t} and
#'         \code{foi_t} must accept arguments \code{(t, par)}; \code{group_pi}
#'         and \code{group_foi} must accept \code{(a, b, par)}.
#'         A warning is issued if \code{group_foi} is provided while \code{t}
#'         is not a two-column matrix.
#'   \item \strong{Allowed values}: \code{foi_functional_form} and
#'         \code{catalytic_model_type} must be among the allowed sets listed
#'         above. If \code{foi_functional_form == "Splines"}, then
#'         \code{catalytic_model_type} must be either \code{"SimpleCatalytic"}
#'         or \code{"WaningImmunity"}.
#'   \item \strong{Initial values}: When any user-supplied modelling functions
#'         are provided, \code{par_init} must be a named numeric vector.
#'   \item \strong{Fixed parameters}: \code{model_fixed_params} must be a list
#'         and is required when using \code{"Griffiths"} (needs \code{tau}),
#'         \code{"PiecewiseConstant"} (needs strictly increasing
#'         \code{upper_cutoffs} whose maximum covers the largest age in \code{t}),
#'         \code{"RestrictedCatalytic"} (needs \code{k}, \code{l} in \[0, 1\]),
#'         and \code{"WaningImmunity"} (needs \code{w} as a nonnegative scalar
#'         or \code{NA} to estimate).
#'   \item \strong{Brent conditions}: When optimisation reduces to one free
#'         parameter (or the constant FOI form under certain presets) and
#'         \code{rho} is provided, Brent requires finite \code{lower} and
#'         \code{upper} with \code{lower < upper}.
#'   \item \strong{rho}: If provided, \code{rho} must lie in \[0, 1\], otherwise
#'         use \code{NA} to request estimation.
#' }
#'
#' On failure the function stops with a single message containing all detected
#' problems. On success it returns \code{TRUE} invisibly.
#'
#' @return \code{TRUE} invisibly on success; otherwise an error is thrown. Some
#'   validations may emit warnings.
#'
#' @seealso \code{\link{FoiFromCatalyticModel}}, \code{\link{set_pi_t}},
#'   \code{\link{set_foi_t}}, \code{\link{set_group_pi}}, \code{\link{set_group_foi}}
#'
#' @keywords internal
#' @noRd
check_inputs <- function(t, y, n, pi_t = NA, foi_t = NA, group_pi = NULL, group_foi = NULL, par_init = NA, rho = 1, catalytic_model_type = NA, foi_functional_form = NA, model_fixed_params = NA, lower = -Inf, upper = Inf) {

  errors  <- character()
  warns   <- character()

  # ---------- helpers ----------
  is_missing_scalar <- function(x) {
    if (is.list(x) || is.function(x)) return(FALSE)
    is.null(x) || (length(x) == 1 && (is.na(x) || identical(x, NA)))
  }
  provided_fun <- function(f) is.function(f)
  provided_any <- function(x) !is_missing_scalar(x)

  has_args <- function(f, needed) {
    if (!is.function(f)) return(FALSE)
    fn <- names(formals(f))
    all(needed %in% fn)
  }

  t_is_vec <- is.atomic(t) && is.numeric(t) && is.null(dim(t))
  t_is_mat <- is.matrix(t) && ncol(t) == 2 && is.numeric(t)
  if (!(t_is_vec || t_is_mat)) {
    errors <- c(errors, "t must be either a numeric vector or a numeric 2-column matrix.")
  }

  # y, n checks
  if (!is.numeric(y) || !is.numeric(n)) {
    errors <- c(errors, "y and n must be numeric vectors.")
  } else {
    expected_len <- if (t_is_mat) nrow(t) else if (t_is_vec) length(t) else NA
    if (!is.na(expected_len)) {
      if (length(y) != expected_len || length(n) != expected_len) {
        errors <- c(errors, sprintf("y and n must each have length %d (to match t).", expected_len))
      }
    }
  }

  # ---------- model choice cross-checks ----------
  # 1) If t is a vector, require exactly one of pi_t OR catalytic_model_type (not both)
  if (t_is_vec) {
    has_pi   <- provided_fun(pi_t)
    has_cat  <- provided_any(catalytic_model_type) && !is_missing_scalar(catalytic_model_type)
    if (has_pi == has_cat) {
      errors <- c(errors, "For vector t: provide exactly ONE of pi_t OR catalytic_model_type (not both).")
    }
  }

  # 2) Must always provide exactly one of foi_t OR foi_functional_form
  has_foi_fun   <- provided_fun(foi_t)
  has_foi_form  <- provided_any(foi_functional_form) && !is_missing_scalar(foi_functional_form)
  if (has_foi_fun == has_foi_form) {
    errors <- c(errors, "Provide exactly ONE of foi_t OR foi_functional_form (not both).")
  }

  # foi_t signature
  if (has_foi_fun && !has_args(foi_t, c("t", "par"))) {
    errors <- c(errors, "foi_t must be a function with arguments (t, par).")
  }

  # Optional group_foi: if provided, must be function (a, b, par)
  if (!is.null(group_foi)) {
    if (!provided_fun(group_foi) || !has_args(group_foi, c("a", "b", "par"))) {
      errors <- c(errors, "group_foi (if provided) must be a function with arguments (a, b, par).")
    }
    if (!t_is_mat) {
      warns <- c(warns, "group_foi is only useful when t is a 2-column matrix of age groups.")
    }
  }

  # 3) Must provide at least one of: pi_t, group_pi, catalytic_model_type
  has_pi_fun   <- provided_fun(pi_t)
  has_gpi_fun  <- !is.null(group_pi) && provided_fun(group_pi)
  has_cat_type <- provided_any(catalytic_model_type) && !is_missing_scalar(catalytic_model_type)

  if (!(has_pi_fun || has_gpi_fun || has_cat_type)) {
    errors <- c(errors, "Provide at least one of: pi_t, group_pi, or catalytic_model_type.")
  }

  # If catalytic_model_type chosen, then pi_t and group_pi should NOT be provided
  if (has_cat_type && (has_pi_fun || has_gpi_fun)) {
    errors <- c(errors, "If catalytic_model_type is provided, do not also provide pi_t or group_pi.")
  }

  # pi_t signature
  if (has_pi_fun && !has_args(pi_t, c("t", "par"))) {
    errors <- c(errors, "pi_t must be a function with arguments (t, par).")
  }

  # group_pi signature
  if (has_gpi_fun && !has_args(group_pi, c("a", "b", "par"))) {
    errors <- c(errors, "group_pi must be a function with arguments (a, b, par).")
  }

  # 4) Allowed sets
  allowed_foi_forms <- c("Constant", "Linear", "Griffiths", "Farringtons", "PiecewiseConstant", "Splines")
  allowed_cat_types <- c("OriginalCatalytic", "RestrictedCatalytic", "SimpleCatalytic", "WaningImmunity", "SimpleCatalytic_NegativeCorrected", "WaningImmunity_NegativeCorrected")

  if (has_foi_form && !(foi_functional_form %in% allowed_foi_forms)) {
    errors <- c(errors, sprintf("foi_functional_form must be one of: %s.",
                                paste(allowed_foi_forms, collapse = ", ")))
  }

  if (has_cat_type && !(catalytic_model_type %in% allowed_cat_types)) {
    errors <- c(errors, sprintf("catalytic_model_type must be one of: %s.",
                                paste(allowed_cat_types, collapse = ", ")))
  }

  # Splines requires catalytic_model_type in {SimpleCatalytic, WaningImmunity}
  if (has_foi_form && foi_functional_form == "Splines") {
    if (!(has_cat_type && catalytic_model_type %in% c("SimpleCatalytic", "WaningImmunity"))) {
      errors <- c(errors, "If foi_functional_form is 'Splines', catalytic_model_type must be 'SimpleCatalytic' or 'WaningImmunity'.")
    }
  }

  # 5) par_init required if any of pi_t, group_pi, foi_t, group_foi are provided
  any_user_fns <- any(c(has_pi_fun, has_gpi_fun, has_foi_fun, !is.null(group_foi) && provided_fun(group_foi)))
  if (any_user_fns) {
    if (is_missing_scalar(par_init) || !is.numeric(par_init) || length(par_init) < 1) {
      errors <- c(errors, "par_init must be provided (numeric vector) when supplying pi_t, group_pi, foi_t, or group_foi.")
    } else {
      if (is.null(names(par_init)) || any(!nzchar(names(par_init)))) {
        errors <- c(errors, "par_init should be a NAMED numeric vector (non-empty names).")
      }
    }
  }

  # 6) model_fixed_params when required
  need_mfp <- FALSE
  if (has_foi_form && foi_functional_form %in% c("Griffiths", "PiecewiseConstant")) need_mfp <- TRUE
  if (has_cat_type && catalytic_model_type %in% c("RestrictedCatalytic", "WaningImmunity")) need_mfp <- TRUE

  if (need_mfp) {
    if (is_missing_scalar(model_fixed_params) || !is.list(model_fixed_params)) {
      errors <- c(errors, "model_fixed_params must be provided as a LIST when required by the chosen model/form.")
    } else {
      # Griffiths: need tau
      if (has_foi_form && foi_functional_form == "Griffiths") {
        if (is.null(model_fixed_params$tau) || !is.numeric(model_fixed_params$tau) || length(model_fixed_params$tau) != 1) {
          errors <- c(errors, "For foi_functional_form='Griffiths', model_fixed_params must contain numeric scalar 'tau' (maternal antibody age cutoff).")
        }
      }
      # PiecewiseConstant: need upper_cutoffs; last >= max age in t
      if (has_foi_form && foi_functional_form == "PiecewiseConstant") {
        uc <- model_fixed_params$upper_cutoffs
        if (is.null(uc) || !is.numeric(uc) || length(uc) < 1) {
          errors <- c(errors, "For 'PiecewiseConstant', model_fixed_params must contain 'upper_cutoffs' (numeric vector).")
        } else {
          if (is.unsorted(uc, strictly = TRUE)) {
            errors <- c(errors, "'upper_cutoffs' must be strictly increasing (e.g., c(10, 20, 60)).")
          }
          max_age <- if (t_is_vec) max(t, na.rm = TRUE) else max(t, na.rm = TRUE)
          if (max(uc) < max_age) {
            errors <- c(errors, sprintf("For 'PiecewiseConstant', max(upper_cutoffs) must be >= largest age in t (%.3f).", max_age))
          }
        }
      }
      # RestrictedCatalytic: need k, l in [0,1]
      if (has_cat_type && catalytic_model_type == "RestrictedCatalytic") {
        k <- model_fixed_params$k; l <- model_fixed_params$l
        if (is.null(k) || is.null(l) || !is.numeric(k) || !is.numeric(l) ||
            length(k) != 1 || length(l) != 1 || k < 0 || k > 1 || l < 0 || l > 1) {
          errors <- c(errors, "For 'RestrictedCatalytic', model_fixed_params must contain k and l in [0, 1].")
        }
      }
      # WaningImmunity: w may be a numeric scalar >= 0, or NA (estimate)
      if (has_cat_type && catalytic_model_type == "WaningImmunity") {
        w <- model_fixed_params$w
        if (is.null(w) || length(w) != 1 || (!is.na(w) && w < 0)) {
          errors <- c(
            errors,
            "For 'WaningImmunity', model_fixed_params$w must be a numeric scalar >= 0, or NA to request estimation."
          )
        }
      }
    }
  }

  # 7) Brent optimisation bounds conditions
  rho_provided <- !is_missing_scalar(rho)
  cond_one_param <- rho_provided && is.numeric(par_init) && length(par_init) == 1
  cond_constant_restricted_simple <- rho_provided && has_foi_form && foi_functional_form == "Constant" &&
    has_cat_type && catalytic_model_type %in% c("RestrictedCatalytic", "SimpleCatalytic")

  if (cond_one_param || cond_constant_restricted_simple) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
      errors <- c(errors, "Using Brent optimisation requires finite 'lower' and 'upper' with lower < upper under the specified conditions.")
    }
  }

  # 8) rho must be NA or in [0,1]
  if (rho_provided && !is.na(rho)) {
    if (!is.numeric(rho) || length(rho) != 1 || rho < 0 || rho > 1) {
      errors <- c(errors, "rho must be NA or a numeric scalar in [0, 1].")
    }
  }

  # ---------- report ----------
  if (length(warns) > 0) warning(paste(warns, collapse = "\n"), call. = FALSE)
  if (length(errors) > 0) stop(paste(c("Input validation failed:", errors), collapse = "\n"), call. = FALSE)

  invisible(TRUE)
}
