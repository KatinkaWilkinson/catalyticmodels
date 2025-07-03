#' Vynnycky Catalytic Model D (Constant FOI, Perfect Assay Sensitivity)
#'
#' Fits a catalytic model assuming a constant force of infection (FOI) across all ages and perfect assay sensitivity (\eqn{ρ = 1}).
#'
#' The seroprevalence function is modeled as:
#' \deqn{
#' \pi(t) = 1 - \exp(-\lambda \cdot (t - 0.5))
#' }
#'
#' This model assumes that all past infections are perfectly detected and estimates a single FOI parameter (\eqn{\lambda}).
#'
#' @param t A numeric vector of age groups.
#' @param y A numeric vector of seropositive counts corresponding to each age group.
#' @param n A numeric vector of total sample sizes corresponding to each age group.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{par}}{Named list with the maximum likelihood estimate for \code{foi} (force of infection).}
#'   \item{\code{CIs}}{Named list with the 95\% bootstrap confidence interval for \code{foi}.}
#'   \item{\code{boot_params}}{Named list containing the bootstrap samples for \code{foi}.}
#' }
#'
#' @export
VynneckyD <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    foi <- par[1]
    pi_t <- 1-exp(-foi*(t-0.5))
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(foi = 0.5)
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 1000
  boot_foi <- numeric(length=boot_num)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n
    bootsamp_params <- tryCatch(
      optim(par = par_init, fn = loglik, t = bootsamp_t, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 1)
    )
    boot_foi[b] <- bootsamp_params[1]
  }
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(foi=params[1]), CIs=list(foi_CI=foi_CI), boot_params=list(boot_foi=boot_foi)))
}
