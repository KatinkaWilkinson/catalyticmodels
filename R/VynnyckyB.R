#' Vynnycky Catalytic Model B (Piecewise FOI, Perfect Assay Sensitivity)
#'
#' Fits a catalytic model with a piecewise constant force of infection (FOI), assuming full assay sensitivity (\eqn{\rho = 1}).
#' This model separates the population into two age groups: younger individuals (\eqn{t < 13}) and older individuals (\eqn{t \geq 13}),
#' each with its own FOI parameter.
#'
#' The seroprevalence function is defined as:
#' \deqn{
#' \pi(t) =
#' \begin{cases}
#'   1 - \exp(-\lambda_y \cdot (t - 0.5)), & \text{for } t < 13 \\
#'   1 - \exp(-\lambda_y \cdot 12.5 - \lambda_o \cdot (t - 13)), & \text{for } t \geq 13
#' \end{cases}
#' }
#'
#' For age intervals \eqn{[a, b]}, the model integrates \eqn{\pi(t)} over the interval and returns the average value.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Maximum likelihood estimates for \code{foi_y} (FOI under age 13) and \code{foi_o} (FOI from age 13 onward).}
#'   \item{CIs}{Bootstrap-based 95% confidence intervals for each parameter.}
#'   \item{boot_params}{Bootstrap samples for \code{foi_y} and \code{foi_o}.}
#' }
#'
#' @details
#' This model captures age-dependent transmission with a structural breakpoint at age 13. It is a special case of Model A,
#' assuming perfect assay sensitivity and that all infections result in detectable seroconversion.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:20
#' y <- round(10 * (1 - exp(-0.3 * pmin(t - 0.5, 12.5) - 0.2 * pmax(0, t - 13))))
#' n <- rep(10, length(t))
#' result <- VynnyckyB(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,2, 2,5, 5,10, 10,15, 15,20), ncol=2, byrow=TRUE)
#' y <- c(1, 3, 6, 8, 9)
#' n <- rep(10, 5)
#' result <- VynnyckyB(age_intervals, y, n)
#' result$CIs
#'
#' @export
VynnyckyB <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    foi_y <- par[1]
    foi_o <- par[2]

    pi_fun <- function(t, foi_y, foi_o) {
      1 - exp(-foi_y * pmin(t - 0.5, 12.5) - foi_o * pmax(0, t - 13))
    }

    if (ncol(t) == 2) {
      a <- t[, 1]
      b <- t[, 2]
      integrate_pi <- function(a, b, foi_y, foi_o) {
        integrand <- function(s) pi_fun(s, foi_y, foi_o)
        avg <- integrate(integrand, lower = a, upper = b)$value / (b - a)
        return(avg)
      }

      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(foi_y = foi_y, foi_o = foi_o))

    } else {
      pi_t <- pi_fun(t, foi_y, foi_o)
    }

    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  par_init <- c(foi_y = 0.5, foi_o = 0.5)
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  boot_num <- 1000
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
      error = function(e) rep(NA, 2)
    )
    boot_foi_y[b] <- bootsamp_params[1]
    boot_foi_o[b] <- bootsamp_params[2]
  }
  foi_y_CI <- quantile(boot_foi_y, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_o_CI <- quantile(boot_foi_o, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(foi_y=params[1], foi_o=params[2]), CIs=list(foi_y_CI=foi_y_CI, foi_o_CI=foi_o_CI), boot_params=list(boot_foi_y=boot_foi_y, boot_foi_o=boot_foi_o)))
}
