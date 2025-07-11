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
#' @param t A one column matrix containing exact values, or a matrix with two columns giving the lower and upper bounds of age intervals.
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
#' t <- matrix(1:10, ncol = 1)
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
#' @importFrom stats optim quantile dbinom
#' @export
MuenchGeneral <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    k <- par[1]
    l <- par[2]
    foi <- par[3]
    if (ncol(t) == 2) { # age buckets
      a <- t[,1]
      b <- t[,2]

      valid <- (b - a) > 0 # if this happens then you will end up dividing by zero below! Therefore check for invalid age buckets
      if (any(!valid)) stop("Invalid interval with zero width detected.")

      pi_t <- k * (l + (exp(-foi*b)-exp(-foi*a)) / (foi*(b-a)) ) # estimate pi_t as the mean of the catalytic function over the range [a,b]. See working out in goodnotes. Comes from 1/(b-a) * integral (a to b) pi(t) dt. I avoided using the integrate function because it is less computationally efficient
    } else { # t represents exact ages (if ncol(t)==1)
      pi_t <- k * (l - exp(-foi*t))
    }
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- dbinom(y, size = n, prob = pi_t, log = TRUE) # I replaced y*log(pi_t) + (n-y)*log(1-pi_t) with dbinom since the function is more robust

    # recall that we are always working with vectors: t, y, n, pi_t, and ll are all vectors. Therefore, the total log likelihood over all the observations = sum(ll)
    return(- sum(ll)) # return negative since optim minimises
  }

  # MLE: maximise the loglik with optim
  par_init <- c(k=0.5, l=1, foi=0.1) # do start values matter? Might be worth looking into, but maybe not the most important thing. Possible solutions: consider adding a check to rerun optim() with different starts if convergence fails or allow the user to optionally provide initial values.
  params <- optim(par = par_init, fn = loglik, t=t, y=y, n=n)$par

  # bootstrap CIs:
  boot_num <- 1000
  boot_k <- numeric(length=boot_num)
  boot_l <- numeric(length=boot_num)
  boot_foi <- numeric(length=boot_num)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)

  # boot_results is a list of length equal to the number of bootstrap samples.
  # Each element is a list of the form list(k, l, foi), returned by optim for that sample.
  # Note: originally I did not check for convergence but then chat pointed this out: Some bootstraps might yield degenerate groups (e.g., all 0 or all 1 in a bin), which makes MLE unstable. Consider adding a convergence check
  boot_results <- lapply(1:boot_num, function(b) {
    boot_samp <- bootstrap_samples[[b]]

    result <- tryCatch(
      optim(par = par_init, fn = loglik, t = boot_samp$t, y = boot_samp$y, n = boot_samp$n),
      error = function(e) NULL
    )

    if (is.null(result) || result$convergence != 0) {
      return(rep(NA, 3))  # return NA if error occurred or convergence failed
    } else {
      return(result$par)  # return the estimated parameters
    }
  })

  # lapply always outputs a list, while sapply (simplify apply, a wrapped function for lapply that simplifies the list output) simplifies the output into a vector/matrix
  boot_k <- sapply(boot_results, function(x) x[[1]]) # pull out all the first elements in each list contained in boot_results, and place all the
  boot_l <- sapply(boot_results, function(x) x[[2]])
  boot_foi <- sapply(boot_results, function(x) x[[3]])

  k_CI <- quantile(boot_k, probs = c(0.025, 0.975), na.rm = TRUE)
  l_CI <- quantile(boot_l, probs = c(0.025, 0.975), na.rm = TRUE)
  foi_CI <- quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(par=list(k=params[1], l=params[2], foi=params[3]), CIs=list(k_CI=k_CI, l_CI=l_CI, foi_CI=foi_CI), boot_params=list(boot_k=boot_k, boot_l=boot_l, boot_foi=boot_foi)))
}
