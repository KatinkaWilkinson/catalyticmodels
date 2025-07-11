#' Muench's Catalytic Model with Constant Force of Infection (Restricted)
#'
#' Fits the classical Muench catalytic model under the assumption of a constant force of infection (FOI),
#' where \eqn{\pi(t) = 1 - \exp(-\lambda t)}. This formulation assumes full susceptibility from birth
#' and that infection accumulates exponentially with age until saturation.
#'
#' For age intervals \eqn{[a, b]}, the model estimates the average seroprevalence within the interval as:
#' \deqn{\bar{\pi}_{[a,b]} = 1 - \frac{e^{-\lambda a} - e^{-\lambda b}}{\lambda(b - a)}}
#' This is derived analytically to avoid computationally expensive numerical integration.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns specifying the lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts corresponding to each age or age interval.
#' @param n A numeric vector of total sample sizes for each age or age interval.
#' @param par_init Optional initial value for \code{foi} (default is 0.1).
#' @param search_range A numeric vector of length 2 giving the lower and upper bounds for optimization via Brent's method (default is c(-1, 3)).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{par}}{Maximum likelihood estimate of the force of infection (\code{foi}).}
#'   \item{\code{CIs}}{Bootstrap-based 95% confidence interval for \code{foi}.}
#'   \item{\code{boot_params}}{Vector of bootstrap estimates for \code{foi}.}
#' }
#'
#' @details
#' This restricted form of the Muench model sets the scaling constants \eqn{k = 1} and \eqn{l = 1}, meaning
#' all individuals are assumed to be fully susceptible and eventually infected if exposed long enough.
#' This makes it appropriate for endemic infections with no significant maternal immunity or heterogeneity in susceptibility.
#' Brent’s method is used for 1D optimization for numerical stability and performance.
#'
#' @examples
#' # Using exact ages
#' t <- matrix(1:10, ncol = 1)
#' y <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' n <- rep(10, 10)
#' result <- MuenchRestricted(t, y, n)
#' result$par
#'
#' # Using age intervals
#' t_int <- matrix(c(0,1, 1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
#' y <- c(0, 1, 2, 3, 4)
#' n <- rep(10, 5)
#' result2 <- MuenchRestricted(t_int, y, n)
#' result2$CIs
#'
#' @importFrom stats optim quantile dbinom
#' @export
MuenchRestricted <- function(t, y, n, par_init = c(foi=0.1), search_range = c(lower = -1, upper = 3)) {
  loglik <- function(par, t, y, n) {
    foi <- par[1]
    if (ncol(t) == 2) { # age buckets
      a <- t[,1]
      b <- t[,2]

      valid <- (b - a) > 0 # if this happens then you will end up dividing by zero below! Therefore check for invalid age buckets
      if (any(!valid)) stop("Invalid interval with zero width detected.")

      pi_t <- 1 + (exp(-foi*b)-exp(-foi*a)) / (foi*(b-a))  # integral of pi_t function over the range [a,b]
    } else { # t represents exact ages (if ncol(t)==1)
      pi_t <- 1 - exp(-foi*t)
    }
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- dbinom(y, size = n, prob = pi_t, log = TRUE) # is this too slow? Should I rather do ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll)) # return negative since optim minimises
  }

  # MLE
  params <- optim(par = par_init, fn = loglik, method = "Brent", lower = search_range[1], upper = search_range[2], t=t, y=y, n=n)$par # according to the optim-generated warning message, Nelder-Mead optimization (the default) is unreliable for 1D optimisation, therefore use Brent for optimising one variable.

  # 95% bootstrap CIs
  boot_num <- 1000
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)

  boot_results <- lapply(1:boot_num, function(b) {
    boot_samp <- bootstrap_samples[[b]]

    result <- tryCatch(
      optim(par = par_init, fn = loglik, t = boot_samp$t, y = boot_samp$y, n = boot_samp$n),
      error = function(e) NULL
    )

    if (is.null(result) || result$convergence != 0) {
      return(rep(NA, 1))  # return NA if error occurred or convergence failed
    } else {
      return(result$par)  # return the estimated parameters
    }
  })

  boot_foi <- sapply(boot_results, function(x) x[[1]])

  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(foi=params[1]), CIs=list(foi_CI=foi_CI), boot_params=list(boot_foi=boot_foi)))
}
