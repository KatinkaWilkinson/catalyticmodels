set_foi_t <- function(type, model_fixed_params = NA) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (type == "MuenchGeneral") {
    foi_t <- function(t, par) {
      foi <- par[3]
      return(foi) # constant foi.
    }
  }

  else if (type == "MuenchRestricted") {
    foi_t <- function(t, par) {
      foi <- par[1]
      return(foi) # constant foi.
    }
  }

  else if (type == "Griffiths") {
    tau <- model_fixed_params$tau
    foi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      return(ifelse(t <= tau, 0, gamma0 * (t + gamma1)))
    }
    return(foi_t)
  }

  else if (type == "Farringtons") {
    foi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      return((gamma0*t + gamma1)*exp(-gamma2*t)+gamma1)
    }
  }

  else if (type == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
    foi_t <- function(t, par) {
      # Each interval is [lower_cutoffs[i], upper_cutoffs[i])
      interval_index <- findInterval(t, lower_cutoffs, rightmost.closed = FALSE)
      return(par[interval_index])
    }
  }

  else if (type == "Splines") {
    foi_t <- function(t, spline_pi_t) {
      pi <- predict(spline_pi_t, t)$y # y is the smooth.spline response name. "y" here is predicting pi(t) (we inputted y/n - this is what we are trying to predict. Remember, y/n=pi(t))
      dpi_t_dt <- predict(spline_pi_t, t, deriv=1)$y
      pi <- pmin(pi, 1 - 1e-8) # to prevent a divide by 0 error
      foi <- dpi_t_dt / (1 - pi)

      # dpi_t_dt <- diff(pi_t)/diff(t)
      # foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)]) this is what I used to have but apparently using splines is better for finding the derivative of pi_t

      return(foi)
    }
  }

  else {
    return(NA)
  }

  return(foi_t)
}

