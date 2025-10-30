#' Fit Catalytic Model and Estimate Force of Infection (sequential)
#'
#' Fits a catalytic model to seroprevalence data to estimate parameters and the
#' force of infection (FOI), with bootstrap confidence intervals. Works for
#' exact-age data or age-interval data and supports preset model choices via
#' \code{catalytic_model_type} and \code{foi_functional_form}, including a spline option.
#' This version runs bootstrap fits sequentially (no parallel backend).
#'
#' @param t Numeric vector of exact ages or a two-column matrix of age intervals \[a, b\].
#' @param y Numeric vector of the number of seropositives in each age or age group.
#' @param n Numeric vector of totals tested in each age or age group.
#' @param pi_t Function. Point-prevalence function of age and parameters. Used when
#'   \code{catalytic_model_type} and \code{foi_functional_form} are not provided.
#' @param foi_t Function. FOI function of age and parameters (or a spline-based FOI when
#'   \code{foi_functional_form == "Splines"}). Used when presets are not provided, and also returned.
#' @param group_pi Function. Group prevalence for an interval \[a, b\] and parameters. Optional; if \code{NULL},
#'   it is set from \code{catalytic_model_type}/\code{foi_functional_form} or from \code{pi_t}.
#' @param group_foi Function. Group FOI for an interval \[a, b\] and parameters. Optional; if \code{NULL},
#'   it is set from \code{catalytic_model_type}/\code{foi_functional_form} or from \code{foi_t}.
#' @param par_init Numeric vector of initial parameter values. If \code{NA}, defaults are set from
#'   \code{catalytic_model_type}/\code{foi_functional_form}.
#' @param rho Numeric scalar or \code{NA}. Test adjustment parameter (for sensitivity/specificity).
#'   If \code{NA}, \code{rho} is estimated with bounds \[0.01, 1\].
#' @param catalytic_model_type Character or \code{NA}. High-level model family selector used to configure
#'   \code{pi_t}, \code{group_pi}, \code{foi_t}, and \code{group_foi} via internal \code{set_*} helpers.
#' @param foi_functional_form Character or \code{NA}. FOI functional form. If \code{"Splines"}, a smoothed
#'   prevalence is fit and FOI is derived from the spline; otherwise preset analytic forms are used.
#' @param model_fixed_params Optional list. Fixed parameter values passed to the internal \code{set_*} helpers.
#' @param boot_num Integer. Number of successful bootstrap fits to retain for confidence intervals (default 1000).
#' @param lower Numeric vector or scalar of lower bounds for parameters (default \code{-Inf}).
#' @param upper Numeric vector or scalar of upper bounds for parameters (default \code{Inf}).
#' @param maxit Integer. Maximum iterations for \code{\link[stats]{optim}} (default 100).
#' @param factr Numeric. Tolerance for \code{L-BFGS-B} (default 1e7). See \code{\link[stats]{optim}}.
#' @param reltol Numeric. Reserved for possible future use in convergence control.
#' @param trace Integer. If > 0, passes a trace flag to \code{optim}.
#' @param convergence_attempts Integer. Maximum number of retries per fit with random jittered starts (default 20).
#' @param tau Numeric. Model-specific hyperparameter passed to \code{set_*} helpers.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item If \code{catalytic_model_type} and/or \code{foi_functional_form} are provided,
#'         it configures \code{pi_t}, \code{group_pi}, \code{foi_t}, and \code{group_foi} via internal helpers.
#'   \item Estimates parameters by maximum likelihood with \code{\link[stats]{optim}}:
#'         \code{"Brent"} for one parameter, otherwise \code{"L-BFGS-B"}. If convergence fails, it retries up to
#'         \code{convergence_attempts} times with jittered initial values clipped to \code{[lower, upper]}.
#'   \item Performs bootstrap resampling of \code{y} conditional on observed \code{t} and \code{n}. A 10% surplus of
#'         bootstrap samples is generated to allow for failed fits; the first \code{boot_num} that converge are kept.
#'   \item Computes FOI by exact age or by age interval. For \code{foi_functional_form == "Splines"},
#'         a smoothed prevalence (\code{\link[stats]{smooth.spline}}) is fit and FOI is derived from it.
#' }
#'
#' @return
#' If \code{foi_functional_form != "Splines"} a list with:
#' \itemize{
#'   \item \code{params_MLE}: Named list of MLEs.
#'   \item \code{params_CI}: Named list of 95% CIs for parameters (percentile).
#'   \item \code{foi_MLE}: Named list of FOI estimates by age or age-interval label.
#'   \item \code{foi_CIs}: Named list of 95% CIs for FOI by age or age-interval label.
#'   \item \code{bootparams}: Matrix of bootstrap parameter estimates (rows = replicates).
#'   \item \code{foi_t}, \code{pi_t}, \code{group_pi}: Functions used, returned for downstream use.
#'   \item \code{tau}: The \code{tau} value provided to the fit.
#' }
#'
#' If \code{foi_functional_form == "Splines"} a list with:
#' \itemize{
#'   \item \code{foi}: Named list of FOI estimates by age or age-interval label.
#'   \item \code{foi_CI}: Named list of 95% CIs for FOI.
#'   \item \code{foi_grid}: FOI evaluated on a regular age grid spanning the data.
#'   \item \code{t}, \code{n}: The input ages (or intervals) and totals.
#'   \item \code{boot_y}: Bootstrap resampled responses used to form the CIs.
#'   \item \code{foi_t}: FOI function used to map from the spline fit to FOI.
#'   \item \code{spline_pi_t}: The fitted \code{\link[stats]{smooth.spline}} object for prevalence.
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item If \code{rho} is \code{NA}, it is estimated with bounds \[0.01, 1\], and \code{lower}/\code{upper} are
#'         extended accordingly.
#'   \item Age-interval labels are formatted as \code{"[a,b)"} in the outputs.
#'   \item This is the sequential version; see \code{FoiFromCatalyticModel()} for the parallel bootstrap variant.
#' }
#'
#' @examples
#' # Age-interval dataset and midpoint exact-age variant
#' t <- matrix(
#'   c(
#'     0, 1,
#'     1, 5,
#'     5, 10,
#'     10, 15,
#'     15, 20,
#'     20, 30,
#'     30, 40,
#'     40, 50,
#'     50, 60
#'   ),
#'   ncol = 2, byrow = TRUE
#' )
#' t_mid <- rowMeans(t)  # c(0.5, 3.0, 7.5, 12.5, 17.5, 25.0, 35.0, 45.0, 55.0)
#' y <- c(28, 224, 577, 712, 740, 808, 848, 873, 894)
#' n <- rep(1000, length(y))
#'
#' # Constant FOI, exact-age midpoints (reduced boot_num for examples)
#' Constant_exactage <- FoiFromCatalyticModel_unparallelised(
#'   t = t_mid, y = y, n = n,
#'   catalytic_model_type = "SimpleCatalytic",
#'   foi_functional_form = "Constant",
#'   lower = 0, upper = 1,
#'   boot_num = 100
#' )
#'
#' # Constant FOI, age-interval data
#' Constant_agegroup <- FoiFromCatalyticModel_unparallelised(
#'   t = t, y = y, n = n,
#'   catalytic_model_type = "SimpleCatalytic",
#'   foi_functional_form = "Constant",
#'   lower = 0, upper = 1,
#'   boot_num = 100
#' )
#'
#' @export
#' @importFrom stats optim quantile smooth.spline runif
FoiFromCatalyticModel_unparallelised <- function(t, y, n, pi_t=NA, foi_t = NA, group_pi = NULL, group_foi = NULL, par_init=NA, rho=1, catalytic_model_type = NA, foi_functional_form = NA, model_fixed_params = NA, boot_num = 1000, lower = -Inf, upper = Inf, maxit = 100, factr = 1e7, reltol = 1e-8, trace = 0, convergence_attempts=20, tau=0) {
  # check that inputs have been entered correctly!
  check_inputs(t, y, n, pi_t, foi_t, group_pi, group_foi,
               par_init, rho, catalytic_model_type,
               foi_functional_form, model_fixed_params,
               lower, upper)

  # Preset: Set pi_t, group_pi, foi_t, group_foi to the correct values, if type != NA
  if (!is.na(foi_functional_form)) {
    foi_t <- set_foi_t(catalytic_model_type, foi_functional_form, model_fixed_params, tau)
    group_foi <- set_group_foi(catalytic_model_type, foi_functional_form, model_fixed_params, foi_t, tau)
  }
  if (!is.na(catalytic_model_type)) {
    if (is.na(foi_functional_form) || !is.na(foi_functional_form) && foi_functional_form != "Splines") {
      pi_t <- set_pi_t(catalytic_model_type, foi_functional_form, model_fixed_params, foi_t=foi_t, tau)
      group_pi <- set_group_pi(catalytic_model_type, foi_functional_form, model_fixed_params, pi_t, tau)
    }
  }

  if (any(is.na(par_init)) || is.na(rho)) { # come back to this...... might be sketchy
    par_init <- set_par_init(catalytic_model_type, foi_functional_form, model_fixed_params, rho, par_init)
    if (is.na(rho)) {
      lower <- c(lower, 0.01)
      upper <- c(upper, 1)
    }
  }

  if (is.null(group_pi)) {
    group_pi <- set_group_pi(catalytic_model_type, foi_functional_form, model_fixed_params, pi_t)
  }

  if (is.null(group_foi)) {
    group_foi <- set_group_foi(catalytic_model_type, foi_functional_form, model_fixed_params, foi_t)
  }

  param_names <- names(par_init)
  if (is.null(param_names)) {
    param_names <- paste0("param", seq_along(par_init))
  }

  # Step 1: Find parameter MLEs
  if (is.na(foi_functional_form) || foi_functional_form != "Splines") {
    convergence_attempt <- 1
    result <- list(convergence = -1)
    par_init_modified <- par_init
    while (result$convergence != 0 && convergence_attempt < convergence_attempts) {
      if (convergence_attempt != 1) {
        rnd <- runif(length(par_init), -0.5, 0.5) * abs(par_init)
        par_init_modified <- par_init + rnd
        par_init_modified <- pmax(pmin(par_init_modified, upper), lower)
        # to make sure that initial vals are still within limits
      }
      if (length(par_init) == 1) {
        result <- optim(par = par_init_modified,
                        fn = neg_total_binom_loglik,
                        method = "Brent",
                        control = list(maxit = maxit),
                        lower = lower,
                        upper = upper,
                        pi_t=pi_t,
                        group_pi=group_pi,
                        group_foi=group_foi,
                        t=t, y=y, n=n,
                        rho=rho, param_names=param_names)
      } else {
        result <- optim(par = par_init_modified,
                        fn = neg_total_binom_loglik,
                        method = "L-BFGS-B",
                        control = list(maxit = maxit, factr = factr),
                        lower = lower,
                        upper = upper,
                        pi_t = pi_t,
                        group_pi = group_pi,
                        group_foi=group_foi,
                        t = t, y = y, n = n,
                        rho = rho, param_names=param_names)
      }
      convergence_attempt <- convergence_attempt + 1
    }

    if (result$convergence != 0) {
      warning("Convergence failed after ", convergence_attempts, " attempts.")
      return(NULL)
    }
    params_MLE <- result$par
    names(params_MLE) <- param_names  # Copy names from par_init ########### !!!!!!!!! does it throw errors if there are no names??
    params_MLE_list <- as.list(params_MLE)

  }

  bootstrap_samples <- create_boot_samps(t, y, n, boot_num * 1.10) # allow 10% to fail if need be

  if (is.na(foi_functional_form) || foi_functional_form != "Splines") {
    boot_results <- list()
    converged_count <- 0
    total_attempts <- 0
    max_attempts <- length(bootstrap_samples)

    i <- 1
    while (converged_count < boot_num && i <= max_attempts) {
      boot_samp <- bootstrap_samples[[i]]

      convergence_attempt <- 1
      result <- list(convergence = -1)
      par_init_modified <- par_init
      while (is.null(result) && convergence_attempt < convergence_attempts || result$convergence != 0 && convergence_attempt < convergence_attempts) {
          if (convergence_attempt != 1) {
            rnd <- runif(length(par_init), -0.5, 0.5) * abs(par_init)
            par_init_modified <- par_init + rnd
            par_init_modified <- pmax(pmin(par_init_modified, upper), lower)
            # to make sure that initial vals are still within limits
          }
          result <- tryCatch({
            if (length(par_init) == 1) {
              optim(par = par_init_modified,
                    fn = neg_total_binom_loglik,
                    method = "Brent",
                    control = list(maxit = maxit, trace=1),
                    lower = lower,
                    upper = upper,
                    pi_t = pi_t,
                    group_pi = group_pi,
                    t = t, y = boot_samp$y, n = n,
                    rho = rho, param_names=param_names)
            } else {
              optim(par = par_init_modified,
                    fn = neg_total_binom_loglik,
                    method = "L-BFGS-B",
                    control = list(maxit = maxit, factr = factr),
                    lower = lower,
                    upper = upper,
                    pi_t = pi_t,
                    group_pi = group_pi,
                    t = t, y = boot_samp$y, n = n,
                    rho = rho, param_names=param_names)
            }
          }, error = function(e) {
            message("Caught error in optim: ", e$message)
            message(boot_samp$y)
            return(NULL)
          })

          convergence_attempt <- convergence_attempt + 1
        }

        total_attempts <- total_attempts + 1
        if (!is.null(result) && result$convergence == 0) {
          boot_results[[converged_count + 1]] <- result$par
          converged_count <- converged_count + 1
          message(paste0("Bootstrap sample ", i, " converged."))
        } else {
          message(paste0("Bootstrap sample ", i, " failed to converge."))

          # if (!is.null(result)) {
          #   message("  Convergence code: ", result$convergence)
          #   message("  Message: ", result$message)
          #   message("  Par: ", paste(round(result$par, 4), collapse = ", "))
          #   message("  Value: ", result$value)
          #   message("  Bootstrap y:", boot_samp$y)
          # }
        }

        i <- i + 1
      }

      if (converged_count < boot_num) {
        warning(paste("Only", converged_count, "bootstrap fits converged out of", boot_num, "required."))
      }

      # Optionally: report how many failed
      cat("Total attempts: ", total_attempts, "\n")
      cat("Failed fits: ", total_attempts - converged_count, "\n")

      if (length(boot_results) == 0) {
        warning("No bootstrap fits converged; cannot compute confidence intervals.")
        return(NULL)
      }
      boot_matrix <- do.call(rbind, boot_results)
      colnames(boot_matrix) <- param_names

      params_CI <- lapply(1:ncol(boot_matrix), function(i) {
        quantile(boot_matrix[, i], probs = c(0.025, 0.975), na.rm = TRUE)
      })

      names(params_CI) <- param_names

    }



  # Step 3: Find FOI MLE for each age group in t
  if (is.null(dim(t)) || ncol(t) == 1) {
    if (!is.na(foi_functional_form) && foi_functional_form == "Splines") {
      spline_pi_t <- smooth.spline(t, y/n)
      foi <- foi_t(t=t, spline_pi_t=spline_pi_t)

      t_grid <- seq(min(t), max(t), length.out = 100)
      foi_grid <- foi_t(t=t_grid, spline_pi_t=spline_pi_t)
    } else {
      foi_MLE <- mapply(foi_t, t, MoreArgs = list(par = params_MLE))
    }
  } else {
    if (!is.na(foi_functional_form) && foi_functional_form == "Splines") {
      t_mid <- (t[,1] + t[,2]) / 2
      spline_pi_t <- smooth.spline(t_mid, y/n)
      foi <- mapply(group_foi, a = t[, 1], b = t[, 2],MoreArgs = list(spline_pi_t=spline_pi_t))

      t_grid <- seq(min(t[,1]), max(t[,2]), length.out = 100)
      foi_grid <- foi_t(t=t_grid, spline_pi_t=spline_pi_t)
    } else {
      foi_MLE <- mapply(group_foi, t[,1], t[,2], MoreArgs = list(par = params_MLE))
    }
  }


  # Step 4: Find FOI bootstrap CI
  # pseudocode:
  # we have boot_matrix. each row contains the par values for one boot sample (therefore 1000 rows)
  # if t is a vector or 1 col matrix, for each value in t, apply foi_t to every row in boot_matrix
  # if t is a 2 col matrix, for each row in t, pull out a and b and apply group_foi to every row in boot_matrix
  # then use quantile to get the 95% CI for the foi of each t category
  if (!is.na(foi_functional_form) && foi_functional_form == "Splines") {
    foi_mat <- t(sapply( # Transpose to make rows = bootstraps, cols = ages
      bootstrap_samples,
      FUN = function(sample) {
        if (is.null(dim(t)) || ncol(t) == 1) {
          boot_spline_pi_t <- smooth.spline(t, sample$y / n)
          foi_t(t = t, spline_pi_t = boot_spline_pi_t)
        } else {
          t_mid <- (t[, 1] + t[, 2]) / 2
          boot_spline_pi_t <- smooth.spline(t_mid, sample$y / n)
          mapply(
            group_foi,
            a = t[, 1],
            b = t[, 2],
            MoreArgs = list(spline_pi_t = boot_spline_pi_t)
          )
        }
      }
    ))
    boot_y <- unlist(bootstrap_samples)
  } else if (is.null(dim(t)) || ncol(t) == 1) {
    foi_mat <- sapply(t, function(age) {
      apply(boot_matrix, 1, foi_t, t = age)
    })
  } else {
    a_vec <- t[, 1]
    b_vec <- t[, 2]
    foi_mat <- mapply(function(a, b) {
      apply(boot_matrix, 1, group_foi, a = a, b = b)
    }, a = a_vec, b = b_vec)
  }

  foi_CIs <- apply(foi_mat, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

  t_labels <- if (is.null(dim(t)) || ncol(t) == 1) {
    as.character(t)  # exact ages
  } else {
    paste0("[", t[,1], ",", t[,2], ")")  # interval labels
  }

  if (!is.na(foi_functional_form) && foi_functional_form == "Splines") {
    foi_list <- as.list(foi)
    names(foi_list) <- t_labels
  } else {
    foi_MLE_list <- as.list(foi_MLE)
    names(foi_MLE_list) <- t_labels
  }

  foi_CI_list <- lapply(seq_len(ncol(foi_CIs)), function(i) foi_CIs[, i])
  names(foi_CI_list) <- t_labels

  #print(boot_matrix)

  if (!is.na(foi_functional_form) && foi_functional_form == "Splines") {
    return(list(foi=foi_list, foi_CI=foi_CI_list, foi_grid=foi_grid, t=t, boot_y=boot_y, n=n, foi_t=foi_t, spline_pi_t=spline_pi_t))
  } # t, boot_y, n are used by R0... is there a better way to do this?


  return(list(params_MLE = params_MLE_list, params_CI=params_CI, foi_MLE=foi_MLE_list, foi_CIs = foi_CI_list, bootparams = boot_matrix, foi_t=foi_t, group_foi=group_foi, pi_t=pi_t, group_pi=group_pi, tau=tau, rho=rho, model_fixed_params=model_fixed_params, t=t, y=y, n=n)) # output foi_t so that it can be used by plot! output pi_t and group_pi for R0!
}
