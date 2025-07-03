#' Vynnycky Catalytic Model A (Piecewise FOI with Assay Sensitivity)
#'
#' Implements a catalytic model with piecewise constant force of infection (FOI) and a scaling factor \eqn{\rho}
#' that accounts for imperfect assay sensitivity or partial immunity. This model captures two distinct transmission phases:
#'
#' \deqn{
#' \pi(t) =
#' \begin{cases}
#'   \rho \left(1 - \exp(-\lambda_y (t - 0.5))\right), & \text{for } t < 13 \\
#'   \rho \left(1 - \exp(-\lambda_y \cdot 12.5 - \lambda_o (t - 13))\right), & \text{for } t \ge 13
#' \end{cases}
#' }
#'
#' For age intervals \eqn{[a, b]}, the model computes the average seroprevalence by numerically integrating \eqn{\pi(t)} across the interval.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with:
#' \describe{
#'   \item{par}{Maximum likelihood estimates for \code{rho}, \code{foi_y} (youth FOI), and \code{foi_o} (older FOI).}
#'   \item{CIs}{Bootstrap-based 95% confidence intervals for each parameter.}
#'   \item{boot_params}{Bootstrap samples for \code{rho}, \code{foi_y}, and \code{foi_o}.}
#' }
#'
#' @details
#' This model is especially useful in seroepidemiological studies where transmission intensity differs by age,
#' and not all infections lead to seroconversion (or detection), which is captured by the \code{rho} parameter.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:20
#' y <- round(10 * (1 - exp(-0.3 * pmin(t - 0.5, 12.5) - 0.2 * pmax(0, t - 13))))
#' n <- rep(10, length(t))
#' result <- VynnyckyA(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,2, 2,5, 5,10, 10,15, 15,20), ncol=2, byrow=TRUE)
#' y <- c(1, 3, 6, 8, 9)
#' n <- rep(10, 5)
#' result <- VynnyckyA(age_intervals, y, n)
#' result$CIs
#'
#' @export
VynnyckyA <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    rho <- par[1]
    foi_y <- par[2]
    foi_o <- par[3]

    pi_fun <- function(t, rho, foi_y, foi_o) {
      rho * (1 - exp(-foi_y * pmin(t - 0.5, 12.5) - foi_o * pmax(0, t - 13)))
    }

    if (is.matrix(t) && ncol(t) == 2) {
      a <- t[, 1]
      b <- t[, 2]
      integrate_pi <- function(a, b, rho, foi_y, foi_o) {
        integrand <- function(s) pi_fun(s, rho, foi_y, foi_o)
        avg <- integrate(integrand, lower = a, upper = b)$value / (b - a)
        return(avg)
      }

      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(rho = rho, foi_y = foi_y, foi_o = foi_o))

    } else {
      pi_t <- pi_fun(t, rho, foi_y, foi_o)
    }

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
