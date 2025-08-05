FoiFromCatalyticModel <- function(t, y, n, pi_t=NA, group_pi = NA, rho=1, w=0, foi_t = NA, group_foi = NA, type = NA, model_fixed_params = NA, boot_num = 1000, par_init=NA, lower = -Inf, upper = Inf, maxit = 100, factr = 1e7, reltol = 1e-8) {
  # Preset: Set pi_t, group_pi, foi_t, group_foi to the correct values, if type != NA
  if (!is.na(type)) {
    if (type != "Splines") {
      pi_t <- set_pi_t(type, model_fixed_params)
      group_pi <- set_group_pi(type, model_fixed_params)
    }
    foi_t <- set_foi_t(type, model_fixed_params)
    group_foi <- set_group_foi(type, model_fixed_params)
    if (is.na(par_init) && type != "Splines") {
      par_init <- set_par_init(type, model_fixed_params)
    }
  }

  # Step 1: Find parameter MLEs
  if (type != "Splines") {
    if (length(par_init) == 1) {
      params_MLE <- optim(par = par_init,
                          fn = neg_total_binom_loglik,
                          method = "Brent",
                          control = list(maxit = maxit),
                          lower = lower,
                          upper = upper,
                          pi_t=pi_t,
                          group_pi=group_pi,
                          t=t, y=y, n=n,
                          rho=rho)$par
    } else {
      params_MLE <- optim(par = par_init,
                          fn = neg_total_binom_loglik,
                          method = "L-BFGS-B", # try gauss seidel #######!!!!!!!!!!!!!!!!!!!!!!!!
                          control = list(maxit = maxit, factr = factr, reltol=reltol),
                          lower = lower,
                          upper = upper,
                          pi_t = pi_t,
                          group_pi = group_pi,
                          t = t, y = y, n = n,
                          rho = rho)$par
    }
    names(params_MLE) <- names(par_init)  # Copy names from par_init
    params_MLE_list <- as.list(params_MLE)

  }

  bootstrap_samples <- create_boot_samps(t, y, n, boot_num * 2)

  if (type != "Splines") {
    boot_results <- list()
    converged_count <- 0
    total_attempts <- 0
    max_attempts <- length(bootstrap_samples)

    i <- 1
    while (converged_count < boot_num && i <= max_attempts) {
      boot_samp <- bootstrap_samples[[i]]

      result <- tryCatch({
        if (length(par_init) == 1) {
          optim(par = par_init,
                fn = neg_total_binom_loglik,
                method = "Brent",
                control = list(maxit = maxit),
                lower = lower,
                upper = upper,
                pi_t = pi_t,
                group_pi = group_pi,
                t = t, y = boot_samp$y, n = n,
                rho = rho)
        } else {
          optim(par = par_init,
                fn = neg_total_binom_loglik,
                method = "L-BFGS-B",
                control = list(maxit = maxit, factr = factr, reltol=reltol, trace=1),
                lower = lower,
                upper = upper,
                pi_t = pi_t,
                group_pi = group_pi,
                t = t, y = boot_samp$y, n = n,
                rho = rho)
        }
      }, error = function(e) {
        message("Caught error in optim: ", e$message)
        return(NULL)
      })

      total_attempts <- total_attempts + 1
      if (!is.null(result) && result$convergence == 0) {
        boot_results[[converged_count + 1]] <- result$par
        converged_count <- converged_count + 1
      } else {
        message(paste0("Bootstrap sample ", i, " failed to converge."))

        if (!is.null(result)) {
          message("  Convergence code: ", result$convergence)
          message("  Message: ", result$message)
          message("  Par: ", paste(round(result$par, 4), collapse = ", "))
          message("  Value: ", result$value)
          message("  Bootstrap y:", boot_samp$y)
        }
      }

      i <- i + 1
    }

    if (converged_count < boot_num) {
      warning(paste("Only", converged_count, "bootstrap fits converged out of", boot_num, "required."))
    }

    # Optionally: report how many failed
    cat("Total attempts: ", total_attempts, "\n")
    cat("Failed fits: ", total_attempts - converged_count, "\n")

    boot_matrix <- do.call(rbind, boot_results)

    params_CI <- lapply(1:ncol(boot_matrix), function(i) {
      quantile(boot_matrix[, i], probs = c(0.025, 0.975), na.rm = TRUE)
    })

    names(params_CI) <- names(par_init)

  }

  # Step 3: Find FOI MLE for each age group in t
  if (is.null(dim(t)) || ncol(t) == 1) {
    if (type == "Splines") {
      spline_pi_t <- smooth.spline(t, y/n)
      foi <- foi_t(t=t, spline_pi_t=spline_pi_t)

      t_grid <- seq(min(t), max(t), length.out = 100)
      foi_grid <- foi_t(t=t_grid, spline_pi_t=spline_pi_t)
    } else {
      foi_MLE <- mapply(foi_t, t, MoreArgs = list(par = params_MLE))
    }
  } else {
    if (type == "Splines") {
      t_mid <- (t[,1] + t[,2]) / 2
      spline_pi_t <- smooth.spline(t_mid, y/n)
      foi <- mapply(group_foi, a = t[, 1], b = t[, 2])

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
  if (type == "Splines") {
    foi_mat <- sapply(bootstrap_samples, function(sample) {
      spline_pi_t <- smooth.spline(t, sample$y / n)
      foi_t(t = t, spline_pi_t = spline_pi_t)
    })
    foi_mat <- t(foi_mat)     # Transpose to make rows = bootstraps, cols = ages
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

  if (type == "Splines") {
    foi_list <- as.list(foi)
    names(foi_list) <- t_labels
  } else {
    foi_MLE_list <- as.list(foi_MLE)
    names(foi_MLE_list) <- t_labels
  }

  foi_CI_list <- lapply(seq_len(ncol(foi_CIs)), function(i) foi_CIs[, i])
  names(foi_CI_list) <- t_labels

  #print(boot_matrix)

  if (type == "Splines") {
    return(list(foi=foi_list, foi_CI=foi_CI_list, foi_grid=foi_grid, boot_y, foi_t, spline_pi_t))
  }

  return(list(params_MLE = params_MLE_list, params_CI=params_CI, foi_MLE=foi_MLE_list, foi_CIs = foi_CI_list, bootparams = boot_matrix, foi_t=foi_t)) # output foi_t so that it can be used by plot!
}


