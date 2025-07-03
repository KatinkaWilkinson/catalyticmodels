#' Vynnycky Catalytic Model C (Constant FOI with Assay Sensitivity)
#'
#' Fits a catalytic model assuming a constant force of infection (FOI) across all ages, while allowing for imperfect assay sensitivity (\eqn{ρ}).
#'
#' The seroprevalence function is modeled as:
#' \deqn{
#' \pi(t) = \rho \cdot \left(1 - \exp(-\lambda \cdot (t - 0.5))\right)
#' }
#'
#' This model allows for a single FOI parameter (\eqn{\lambda}) and a sensitivity parameter (\eqn{\rho}) representing the proportion of past infections correctly detected by the assay.
#'
#' @param t A numeric vector of age groups.
#' @param y A numeric vector of seropositive counts corresponding to each age group.
#' @param n A numeric vector of total sample sizes corresponding to each age group.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{par}}{Named list of maximum likelihood estimates: \code{rho} (assay sensitivity) and \code{foi} (force of infection).}
#'   \item{\code{CIs}}{Named list of 95\% bootstrap confidence intervals for \code{rho} and \code{foi}.}
#'   \item{\code{boot_params}}{Named list containing the bootstrap samples for \code{rho} and \code{foi}.}
#' }
#'
#' @export
VynnyckyC <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    rho <- par[1]
    foi <- par[2]
    pi_t <- rho*(1-exp(-foi*(t-0.5)))
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(rho=0.8,foi = 0.5)
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 1000
  boot_rho <- numeric(length=boot_num)
  boot_foi <- numeric(length=boot_num)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n
    bootsamp_params <- tryCatch(
      optim(par = par_init, fn = loglik, t = bootsamp_t, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 2)
    )
    boot_rho[b] <- bootsamp_params[1]
    boot_foi[b] <- bootsamp_params[2]
  }
  rho_CI <- quantile(boot_rho, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(rho=params[1], foi=params[2]), CIs=list(rho_CI=rho_CI, foi_CI=foi_CI), boot_params=list(boot_rho=boot_rho, boot_foi=boot_foi)))
}
