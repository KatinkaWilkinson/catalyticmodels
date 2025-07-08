#' Muench's Catalytic Model with Constant Force of Infection (Restricted)
#'
#' Fits the classical Muench catalytic model under the assumption of a constant force of infection (FOI),
#' where \eqn{\pi(t) = 1 - \exp(-\lambda t)}. This form assumes full susceptibility from birth,
#' and that infection risk accumulates exponentially with age.
#'
#' For age intervals \eqn{[a, b]}, the model computes the average of \eqn{\pi(t)} over the interval:
#' \deqn{\bar{\pi}_{[a,b]} = \frac{1}{b - a} \int_a^b \left(1 - e^{-\lambda t}\right) dt = 1 - \frac{e^{-\lambda a} - e^{-\lambda b}}{\lambda(b - a)}}
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{par}: Maximum likelihood estimate for the force of infection \code{foi}.
#'   \item \code{CIs}: Bootstrap-based 95% confidence interval for \code{foi}.
#'   \item \code{boot_params}: Bootstrap samples for \code{foi}.
#' }
#'
#' @details
#' This restricted form of the catalytic model fixes the scaling constants \eqn{k = 1} and \eqn{l = 1},
#' implying that all individuals are equally susceptible and that infection eventually reaches saturation.
#' It is suitable for pathogens with long-term endemic circulation and no maternal immunity window.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:10
#' y <- c(0, 1, 1, 2, 3, 4, 5, 6, 7, 8)
#' n <- rep(10, length(t))
#' result <- MuenchRestricted(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,1, 1,2, 2,3, 3,4, 4,5), ncol=2, byrow=TRUE)
#' y <- c(0, 1, 2, 3, 4)
#' n <- rep(10, 5)
#' result <- MuenchRestricted(age_intervals, y, n)
#'
#' @export
MuenchRestricted <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    foi <- par[1]
    if (ncol(t) == 2) { # age buckets
      a <- t[,1]
      b <- t[,2]
      pi_t <- 1 + (exp(-foi*b)-exp(-foi*a)) / (foi*(b-a))  # integral of pi_t function over the range [a,b]
    } else { # t is a vector
      pi_t <- 1 - exp(-foi*t)
    }
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll)) # return negative since optim minimises
  }

  # MLE
  par_init <- c(foi=0.1)
  params <- optim(par = par_init, fn = loglik, method = "Brent", lower = 0, upper = 5, t=t, y=y, n=n)$par # according to the optim-generated warning message, Nelder-Mead optimization (the default) is unreliable for 1D optimisation, therefore use Brent for optimising one variable. - DID I CHOOSE THE LOWER AND UPPER BOUNDS CORRECTLY>???

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
      optim(par = par_init, fn = loglik, method = "Brent", lower = 0, upper = 5, t = bootsamp_t, y = bootsamp_y, n = bootsamp_n)$par,
      error = function(e) rep(NA, 1)
    )
    boot_foi[b] <- bootsamp_params[1]
  }
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(foi=params[1]), CIs=list(foi_CI=foi_CI), boot_params=list(boot_foi=boot_foi)))
}
