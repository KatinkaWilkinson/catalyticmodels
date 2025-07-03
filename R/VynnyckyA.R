#' Vynnycky Catalytic Model A
#'
#' Fits a piecewise catalytic model of seroprevalence that accounts for different forces of infection
#' in younger (<13 years) and older (≥13 years) individuals, and includes an unknown assay sensitivity parameter.
#'
#' The model assumes:
#' \deqn{π(t) = ρ(1 - e^{-λ_y (t - 0.5)})} for \eqn{t < 13}
#' \deqn{π(t) = ρ(1 - e^{-λ_y (12.5) - λ_o (t - 13)})} for \eqn{t \geq 13}
#' where:
#' \itemize{
#'   \item \eqn{ρ} is the assay sensitivity (maximum achievable seroprevalence),
#'   \item \eqn{λ_y} is the force of infection for individuals under 13,
#'   \item \eqn{λ_o} is the force of infection for individuals 13 and older.
#' }
#'
#' @param t A numeric vector of age groups.
#' @param y A numeric vector of seropositive counts corresponding to each age group.
#' @param n A numeric vector of total sample sizes corresponding to each age group.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{par}}{Named list of MLE estimates for \code{rho}, \code{foi_y}, and \code{foi_o}.}
#'   \item{\code{CIs}}{Named list of 95% bootstrap confidence intervals for each parameter.}
#'   \item{\code{boot_params}}{Named list of bootstrap samples for each parameter. Useful for diagnostics or visualisation.}
#' }
#'
#' @export
VynnyckyA <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    rho <- par[1]
    foi_y <- par[2]
    foi_o <- par[3]
    pi_t <- rho*(1-exp(-foi_y*pmin(t-0.5, 12.5) -foi_o*pmax(0,t-13)))
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(rho=0.8,foi_y = 0.5, foi_o = 0.5)
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 1000
  boot_rho <- numeric(length=boot_num)
  boot_foi_y <- numeric(length=boot_num)
  boot_foi_o <- numeric(length=boot_num)
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
    boot_rho[b] <- bootsamp_params[1]
    boot_foi_y[b] <- bootsamp_params[2]
    boot_foi_o[b] <- bootsamp_params[3]
  }
  rho_CI <- quantile(boot_rho, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_y_CI <- quantile(boot_foi_y, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_o_CI <- quantile(boot_foi_o, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(rho=params[1], foi_y=params[2], foi_o=params[3]), CIs=list(rho_CI=rho_CI, foi_y_CI=foi_y_CI, foi_o_CI=foi_o_CI), boot_params=list(boot_rho=boot_rho, boot_foi_y=boot_foi_y, boot_foi_o=boot_foi_o)))
}
