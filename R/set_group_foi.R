set_group_foi <- function(catalytic_model_type, foi_functional_form, model_fixed_params = NA, foi_t = NULL) { # type is a string, model_fixed_params is a list
  # Handles MuenchGeneral, MuenchRestricted, Griffiths, Farringtons, PiecewiseConstant, Splines, Keidings
  if (!is.na(foi_functional_form)) {
    if (foi_functional_form == "Constant") {
      group_foi <- function(a, b, par) {
        foi <- par[["foi"]]
        return(foi) # constant foi
      }
    }

    if (foi_functional_form == "Linear") {
      group_foi <- function(a, b, par) {
        m <- par[["m"]]
        c <- par[["c"]]

        return((m*a + c + m*b + c)/2)
      }
    }


    else if (foi_functional_form == "Griffiths") {
      tau <- model_fixed_params$tau
      group_foi <- function(a, b, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]

        if (b <= tau) return(0)

        # Adjust bounds if partially below tau
        a_adj <- max(a, tau)

        avg_foi<-gamma0*((b-a)/2 + a) + gamma0*gamma1
        return(avg_foi)
      }
    }

    else if (foi_functional_form == "Farringtons") {
      group_foi <- function(a, b, par) {
        gamma0 <- par[["gamma0"]]
        gamma1 <- par[["gamma1"]]
        gamma2 <- par[["gamma2"]]
        if (abs(gamma2) < 1e-6) gamma2 <- 1e-6
        integral <- (((gamma1 - b * gamma0) * gamma2 - gamma0) * exp(-b * gamma2) + exp(-a * gamma2) * ((b - a) * gamma1 * gamma2^2 * exp(a * gamma2) + (a * gamma0 - gamma1) * gamma2 + gamma0)) / gamma2^2
        avg_foi <- integral / (b - a)
        return(avg_foi)
      }
    }

    else if (foi_functional_form == "PiecewiseConstant") {
      upper_cutoffs <- model_fixed_params$upper_cutoffs
      lower_cutoffs <- c(0, upper_cutoffs[-length(upper_cutoffs)])
      group_foi <- function(a, b, par) {
        par <- par[ names(par) != "rho" ]
        # Determine overlap between each piecewise interval and [a, b)
        overlaps <- pmin(b, upper_cutoffs) - pmax(a, lower_cutoffs)
        overlaps[overlaps < 0] <- 0  # No overlap = set to zero
        total_length <- b - a
        weighted_avg <- sum(par * overlaps) / total_length
        return(weighted_avg)
      }
    }

    else if (!is.na(catalytic_model_type) && catalytic_model_type == "SimpleCatalytic" && foi_functional_form == "Splines") {
      group_foi <- function(a, b, spline_pi_t) {
        integrand <- function(x) {
          pi_x <- predict(spline_pi_t, x)$y
          dpi_x <- predict(spline_pi_t, x, deriv = 1)$y
          pi_x <- pmin(pi_x, 1 - 1e-8) # to prevent a divide by 0 error
          dpi_x / (1 - pi_x) # estimated foi
        }
        result <- tryCatch({
          integrate(integrand, a, b)$value / (b - a) #  approximate the mean FOI over [a, b)  by integrating the FOI (π′ / (1−π)) over that interval and dividing by the interval width
        }, error = function(e) NA)
        return(result)
      }
    }

    else if (!is.na(catalytic_model_type) && catalytic_model_type == "WaningImmunity" && foi_functional_form == "Splines") {
      w <- model_fixed_params$w
      group_foi <- function(a, b, spline_pi_t) {
        integrand <- function(x) {
          pi_x <- predict(spline_pi_t, x)$y
          dpi_x <- predict(spline_pi_t, x, deriv = 1)$y
          pi_x <- pmin(pi_x, 1 - 1e-8) # to prevent a divide by 0 error
          (dpi_x + pi_x*w) / (1 - pi_x)
        }
        result <- tryCatch({
          integrate(integrand, a, b)$value / (b - a) #  approximate the mean FOI over [a, b)  by integrating the FOI (π′ / (1−π)) over that interval and dividing by the interval width
        }, error = function(e) NA)
        return(result)
      }
    }

  }

  else {
    group_foi <- function(a, b, par) {
      1/(b-a) * integrate(function(x) {foi_t(x, par)}, a, b)$value
    }
  }

  return(group_foi)
}

