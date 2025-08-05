set_pi_t <- function(type, model_fixed_params = NA) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (type == "MuenchGeneral") {
    pi_t <- function(t, par) {
      k <- par[1]
      l <- par[2]
      foi <- par[3]
      return(k * (l - exp(-foi * t)))
    }
  }

  else if (type == "MuenchRestricted") {
    k <- model_fixed_params$k
    l <- model_fixed_params$l
    pi_t <- function(t, par) {
      foi <- par[1]
      return(k * (l - exp(-foi * t)))
    }
  }

  else if (type == "Griffiths") {
    tau <- model_fixed_params$tau
    pi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      return(ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau)))))
    }
  }

  else if (type == "Farringtons") {
    pi_t <- function(t, par) {
      gamma0 <- par[1]
      gamma1 <- par[2]
      gamma2 <- par[3]
      return(1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2)))
    }
  }

  else if (type == "PiecewiseConstant") {
    upper_cutoffs <- model_fixed_params$upper_cutoffs
    lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])

    pi_t <- function(t, par) {
      foi_pieces <- par
      interval_lengths <- pmax(pmin(t, upper_cutoffs) - lower_cutoffs, 0)
      cum_foi <- sum(foi_pieces * interval_lengths)
      return(1 - exp(-cum_foi))
    }
  }

  else {
    return(NA)
  }

  return(pi_t)
}
