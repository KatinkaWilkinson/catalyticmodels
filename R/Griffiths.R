#' Griffiths' Age-Structured FOI Model with Maternal Antibody Cutoff τ
#'
#' Fits a catalytic model in which the force of infection (FOI) becomes active only after a cutoff age \eqn{\tau}, representing
#' the age beyond which maternal antibodies are no longer protective. The FOI increases linearly with age after this cutoff.
#'
#' The FOI is defined as: \eqn{\lambda(t) = \gamma_0(t + \gamma_1)\mathbf{1}\{t > \tau\}}.
#' The cumulative probability of seroconversion is then:
#' \deqn{\pi(t) = 1 - \exp\left(-\int_0^t \lambda(w) dw\right)}
#'
#' For exact ages, the model directly evaluates \eqn{\pi(t)}. For age intervals \[a, b\], the model integrates \eqn{\pi(t)} over the interval
#' and returns the average value.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#' @param tau The age (in the same units as \code{t}) after which maternal antibodies wane and individuals become susceptible to infection.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Maximum likelihood estimates for \code{gamma0} and \code{gamma1}.}
#'   \item{CIs}{Bootstrap-based 95% confidence intervals for \code{gamma0} and \code{gamma1}.}
#'   \item{boot_params}{Bootstrap samples for \code{gamma0} and \code{gamma1}.}
#' }
#'
#' @details
#' This model is suitable when infection dynamics are known to be negligible during early life due to maternal immunity,
#' and begin increasing with age thereafter.
#'
#' @importFrom stats optim quantile integrate dbinom
#' @export
Griffiths <- function(t, y, n, tau) {
  loglik <- function(par, tau, t, y, n) {
    gamma0 <- par[1]
    gamma1 <- par[2]
    if (ncol(t) == 2) { # in this case we are working with a matrix
      a <- t[,1]
      b <- t[,2]

      pi_fun <- function(t, gamma0, gamma1, tau) {
        out <- ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau))))
        return(out)
      }

      integrate_pi <- function(a, b, gamma0, gamma1, tau) {
        if (b <= tau) {
          return(0)
        } else if (a < tau && b > tau) {
          integrate(pi_fun, lower = tau, upper = b, gamma0 = gamma0, gamma1 = gamma1, tau = tau)$value
        } else {
          integrate(pi_fun, lower = a, upper = b, gamma0 = gamma0, gamma1 = gamma1, tau = tau)$value
        }
      }

      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(gamma0 = gamma0, gamma1 = gamma1, tau = tau)) / (b - a)
    } else { # t is a vector
      pi_t <- pi_t <- 1 - exp((-gamma0*(t^2-tau^2)/2 - gamma0*gamma1*(t-tau))*(t>tau))
    }

    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(gamma0 = 0.1, gamma1 = 1)
  params <- optim(par=par_init, fn=loglik, tau=tau, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 100
  boot_gamma0 <- numeric(length=boot_num)
  boot_gamma1 <- numeric(length=boot_num)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n
    bootsamp_params <- tryCatch(
      optim(par = par_init, fn = loglik, t = bootsamp_t, tau=tau, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 2)
    )
    boot_gamma0[b] <- bootsamp_params[1]
    boot_gamma1[b] <- bootsamp_params[2]
  }
  gamma0_CI <- quantile(boot_gamma0, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma1_CI <- quantile(boot_gamma1, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(gamma0=params[1], gamma1=params[2]), CIs=list(gamma0_CI=gamma0_CI, gamma1_CI=gamma1_CI), boot_params=list(boot_gamma0=boot_gamma0, boot_gamma1=boot_gamma1)))
}
