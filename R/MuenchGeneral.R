#' General Muench Catalytic Model with Scaling Parameters
#'
#' Fits the general form of Muench's catalytic model for seroprevalence, allowing for flexible scaling parameters \eqn{k} and \eqn{l}.
#' The model assumes the cumulative probability of infection by age \eqn{t} is:
#' \deqn{\pi(t) = k(l - \exp(-\lambda t))}
#'
#' This form allows for partial susceptibility (\eqn{k < 1}) and incomplete saturation (\eqn{l < 1}) of infection, making it suitable for modeling diseases with non-universal exposure or immunity.
#'
#' For age intervals \eqn{[a, b]}, the model computes the average value of \eqn{\pi(t)} over the interval using:
#' \deqn{\bar{\pi}_{[a,b]} = k \left(l - \frac{e^{-\lambda a} - e^{-\lambda b}}{\lambda(b - a)}\right)}
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns giving the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{par}: Maximum likelihood estimates for \code{k}, \code{l}, and the force of infection \code{foi}.
#'   \item \code{CIs}: Bootstrap-based 95% confidence intervals for each parameter.
#'   \item \code{boot_params}: Bootstrap samples for \code{k}, \code{l}, and \code{foi}.
#' }
#'
#' @details
#' This general form of the catalytic model relaxes assumptions of complete susceptibility and perfect saturation.
#' It is useful in real-world settings where heterogeneity in exposure or immunity may reduce the total proportion infected over time.
#'
#' @examples
#' # Example with exact ages
#' t <- 1:10
#' y <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' n <- rep(10, length(t))
#' result <- MuenchGeneral(t, y, n)
#' result$par
#'
#' # Example with age intervals
#' age_intervals <- matrix(c(0,1, 1,2, 2,3, 3,4, 4,5), ncol=2, byrow=TRUE)
#' y <- c(0, 1, 2, 3, 4)
#' n <- rep(10, 5)
#' result <- MuenchGeneral(age_intervals, y, n)
#'
#' @export
MuenchGeneral <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    k <- par[1]
    l <- par[2]
    foi <- par[3]
    if (is.matrix(t) && ncol(t) == 2) { # in this case we are working with a matrix
      a <- t[,1]
      b <- t[,2]
      pi_t <- k * (l + (exp(-foi*b)-exp(-foi*a)) / (foi*(b-a)) ) # integral of pi_t function over the range [a,b]
    } else { # t is a vector
      pi_t <- k * (l - exp(-foi*t))
    }
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- dbinom(y, size = n, prob = pi_t, log = TRUE) # I replaced y*log(pi_t) + (n-y)*log(1-pi_t) with dbinom since the function is more robust

    return(- sum(ll)) # return negative since optim minimises
  }

  # MLE: maximise the loglik with optim
  par_init <- c(k=0.5, l=1, foi=0.1)
  params <- optim(par = par_init, fn = loglik, t=t, y=y, n=n)$par

  # bootstrap CIs:
  boot_num <- 1000
  boot_k <- numeric(length=boot_num)
  boot_l <- numeric(length=boot_num)
  boot_foi <- numeric(length=boot_num)
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
    boot_k[b] <- bootsamp_params[1]
    boot_l[b] <- bootsamp_params[2]
    boot_foi[b] <- bootsamp_params[3]
  }
  k_CI <- quantile(boot_k, probs = c(0.025, 0.975), na.rm = TRUE)
  l_CI <- quantile(boot_l, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(k=params[1], l=params[2], foi=params[3]), CIs=list(k_CI=k_CI, l_CI=l_CI, foi_CI=foi_CI), boot_params=list(boot_k=boot_k, boot_l=boot_l, boot_foi=boot_foi)))
}
