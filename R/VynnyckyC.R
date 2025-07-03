#' Vynnycky Catalytic Model C (Constant FOI with Assay Sensitivity)
#'
#' Implements a catalytic model assuming a constant force of infection (FOI) across all ages,
#' while allowing for imperfect assay sensitivity \eqn{\rho}.
#'
#' The seroprevalence function is defined as:
#' \deqn{
#' \pi(t) = \rho \left(1 - \exp(-\lambda (t - 0.5))\right)
#' }
#'
#' This structure includes a sensitivity parameter (\eqn{\rho}) that accounts for imperfect detection or partial immunity,
#' and a constant FOI parameter (\eqn{\lambda}) representing the infection rate across all ages.
#'
#' For age intervals \eqn{[a, b]}, the model computes the average seroprevalence by numerically integrating \eqn{\pi(t)} over the interval.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Maximum likelihood estimates for \code{rho} (assay sensitivity) and \code{foi} (force of infection).}
#'   \item{CIs}{Bootstrap-based 95% confidence intervals for each parameter.}
#'   \item{boot_params}{Bootstrap samples for \code{rho} and \code{foi}.}
#' }
#'
#' @details
#' This model is suitable in settings where the transmission intensity is assumed to be stable across ages
#' and not all infections result in detectable antibody responses.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:15
#' y <- round(10 * (1 - exp(-0.4 * (t - 0.5))))
#' n <- rep(10, length(t))
#' result <- VynnyckyC(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,2, 2,5, 5,10, 10,15), ncol=2, byrow=TRUE)
#' y <- c(1, 3, 6, 8)
#' n <- rep(10, 4)
#' result <- VynnyckyC(age_intervals, y, n)
#' result$CIs
#'
#' @export
VynnyckyC <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    rho <- par[1]
    foi <- par[2]

    pi_fun <- function(t, rho, foi) {
      rho * (1 - exp(-foi * (t - 0.5)))
    }

    if (is.matrix(t) && ncol(t) == 2) {
      a <- t[, 1]
      b <- t[, 2]
      integrate_pi <- function(a, b, rho, foi) {
        integrand <- function(s) pi_fun(s, rho, foi)
        avg <- integrate(integrand, lower = a, upper = b)$value / (b - a)
        return(avg)
      }

      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(rho = rho, foi = foi))

    } else {
      pi_t <- pi_fun(t, rho, foi)
    }

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
