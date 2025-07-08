#' Vynnycky Catalytic Model D (Constant FOI, Perfect Assay Sensitivity)
#'
#' Implements a catalytic model that assumes a constant force of infection (FOI) across all ages and
#' perfect assay sensitivity (i.e., all past infections are detected, with \eqn{ρ = 1}).
#'
#' The seroprevalence function is:
#' \deqn{
#' \pi(t) = 1 - \exp(-\lambda (t - 0.5))
#' }
#'
#' For age intervals \eqn{[a, b]}, the model integrates \eqn{\pi(t)} across the interval and returns the average value.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Maximum likelihood estimate for \code{foi} (force of infection).}
#'   \item{CIs}{Bootstrap-based 95% confidence interval for \code{foi}.}
#'   \item{boot_params}{Bootstrap samples for \code{foi}.}
#' }
#'
#' @details
#' This model is a simplified version of the Muench catalytic model, suitable when assay sensitivity is known to be perfect.
#' It is appropriate in settings where all infections are detectable and transmission is age-independent.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:15
#' y <- round(10 * (1 - exp(-0.3 * (t - 0.5))))
#' n <- rep(10, length(t))
#' result <- VynnyckyD(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,2, 2,5, 5,10, 10,15), ncol=2, byrow=TRUE)
#' y <- c(1, 3, 6, 8)
#' n <- rep(10, 4)
#' result <- VynnyckyD(age_intervals, y, n)
#' result$CIs
#'
#' @export
VynnyckyD <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    foi <- par[1]

    pi_fun <- function(t, foi) {
      1 - exp(-foi * (t - 0.5))
    }

    if (ncol(t) == 2) {
      a <- t[, 1]
      b <- t[, 2]
      integrate_pi <- function(a, b, foi) {
        integrand <- function(s) pi_fun(s, foi)
        avg <- integrate(integrand, lower = a, upper = b)$value / (b - a)
        return(avg)
      }

      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(foi = foi))

    } else {
      pi_t <- pi_fun(t, foi)
    }

    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(foi = 0.5)
  params <- optim(par=par_init, fn=loglik, method = "Brent", lower = -1, upper = 3, t=t, y=y, n=n)$par

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
      optim(par = par_init, fn = loglik, method = "Brent", lower = -1, upper = 3, t = bootsamp_t, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 1)
    )
    boot_foi[b] <- bootsamp_params[1]
  }
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(foi=params[1]), CIs=list(foi_CI=foi_CI), boot_params=list(boot_foi=boot_foi)))
}
