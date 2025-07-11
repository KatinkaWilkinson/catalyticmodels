#' Farrington’s Catalytic Model for Seroprevalence
#'
#' Fits a flexible catalytic model in which the force of infection (FOI) is governed by three parameters: \eqn{\gamma_0}, \eqn{\gamma_1}, and \eqn{\gamma_2}.
#' This model allows the FOI to vary non-linearly with age and includes saturation behavior.
#'
#' The force of infection function is:
#' \deqn{\lambda(t) = (\gamma_0 t + \gamma_1) \exp(-\gamma_2 t) + \gamma_1}
#'
#' The cumulative probability of infection by age \eqn{t} is:
#' \deqn{\pi(t) = 1 - \exp\left(-\int_0^t \lambda(s) ds\right)}
#'
#' This model is suitable when the FOI is expected to increase initially and then decline, reflecting waning exposure or susceptibility with age.
#'
#' @param t A numeric vector of exact ages or midpoints of age intervals.
#' @param y A numeric vector of seropositive counts for each age or age group.
#' @param n A numeric vector of total sample sizes for each age or age group.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{par}}{Maximum likelihood estimates of \code{gamma0}, \code{gamma1}, and \code{gamma2}.}
#'   \item{\code{CIs}}{Bootstrap-based 95% confidence intervals for each parameter.}
#'   \item{\code{boot_params}}{Bootstrap samples for each parameter.}
#' }
#'
#' @details
#' The model estimates the seroprevalence curve by integrating a non-linear FOI function over age. Confidence intervals are computed via nonparametric bootstrap resampling using the `create_boot_samps()` function.
#'
#' If the input `t` is based on age intervals, ensure midpoints are used or modify the function to integrate over intervals.
#'
#' @importFrom stats optim quantile integrate dbinom
#' @export
Farringtons <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    gamma0 <- par[1]
    gamma1 <- par[2]
    gamma2 <- par[3]
    pi_t <- sapply(t, function(ti) {
      integrand <- function(s) (gamma0*s + gamma1)*exp(-gamma2*s)+gamma1
      val <- integrate(integrand, 0, ti)$value
      1 - exp(-val)
    })
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(gamma0 = 0.1, gamma1 = 1, gamma2 = 0.1)
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 1000
  boot_gamma0 <- numeric(length=boot_num)
  boot_gamma1 <- numeric(length=boot_num)
  boot_gamma2 <- numeric(length=boot_num)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n
    bootsamp_params <- tryCatch(
      optim(par = par_init, fn = loglik, t = bootsamp_t, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 3)
    )
    boot_gamma0[b] <- bootsamp_params[1]
    boot_gamma1[b] <- bootsamp_params[2]
    boot_gamma2[b] <- bootsamp_params[3]
  }
  gamma0_CI <- quantile(boot_gamma0, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma1_CI <- quantile(boot_gamma1, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma2_CI <- quantile(boot_gamma2, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(gamma0=params[1], gamma1=params[2], gamma2=params[3]), CIs=list(gamma0_CI=gamma0_CI, gamma1_CI=gamma1_CI, gamma2_CI=gamma2_CI), boot_params=list(boot_gamma0=boot_gamma0, boot_gamma1=boot_gamma1, boot_gamma2=boot_gamma2)))
}
