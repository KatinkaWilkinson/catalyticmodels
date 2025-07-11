#' Griffiths' Age-Structured FOI Model with Maternal Antibody Cutoff τ
#'
#' Fits a catalytic model where the force of infection (FOI) becomes active only after a cutoff age \eqn{\tau},
#' representing the waning of maternal antibodies. The FOI increases linearly with age after this cutoff.
#'
#' The FOI is defined as:
#' \deqn{\lambda(t) = \gamma_0(t + \gamma_1)\mathbf{1}\{t > \tau\}}
#'
#' The cumulative probability of seroconversion up to age \eqn{t} is:
#' \deqn{\pi(t) = 1 - \exp\left(-\int_0^t \lambda(w) dw\right)}
#'
#' For exact ages, the model directly evaluates \eqn{\pi(t)}. For age intervals from \code{a} to \code{b}, it computes the average
#' \eqn{\pi(t)} over the interval, which is used as the binomial success probability.
#'
#' @param t A matrix of ages. For exact ages, use a one-column matrix; for age intervals, use a matrix with two columns giving lower and upper bounds.
#' @param y A numeric vector of seropositive counts corresponding to each age group or interval.
#' @param n A numeric vector of total individuals sampled in each age group or interval.
#' @param tau A positive numeric value specifying the age (in the same units as \code{t}) at which maternal antibodies wane.
#' @param par_init Optional named numeric vector of initial values for the parameters \code{gamma0} and \code{gamma1}.
#' @param boot_num Number of bootstrap replicates to use for computing confidence intervals. Default is 1000.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{par}}{Maximum likelihood estimates for \code{gamma0}, \code{gamma1}, and estimated FOI values for each age or interval.}
#'   \item{\code{CIs}}{Bootstrap-based 95% confidence intervals for \code{gamma0}, \code{gamma1}, and the FOI values per age or interval.}
#'   \item{\code{boot_params}}{Bootstrap samples of \code{gamma0} and \code{gamma1}.}
#' }
#'
#' @details
#' This model is appropriate in settings where maternal immunity provides early protection, and susceptibility increases
#' linearly with age once this protection wanes. FOI estimates and their confidence intervals are returned per age category,
#' based on either exact ages or age intervals. The binomial log-likelihood is maximized, and uncertainty is estimated using a nonparametric bootstrap.
#'
#' The number of bootstrap replicates used to compute parameter and FOI confidence intervals can be controlled with the \code{boot_num} argument.
#'
#' @note This function requires a helper function \code{create_boot_samps()} that returns bootstrap samples in the form:
#' \code{list(list(t = ..., y = ..., n = ...), ...)}.
#'
#' @examples
#' # Example with exact ages
#' t <- matrix(c(1, 3, 5, 7, 9), ncol = 1)
#' y <- c(0, 3, 6, 9, 10)
#' n <- c(10, 10, 10, 10, 10)
#' tau <- 2
#' result <- Griffiths(t, y, n, tau, boot_num = 30) # kept low so that example can run quickly
#' result$par
#' result$CIs
#'
#' # Example with age intervals
#' t_int <- matrix(c(0,2, 2,5, 5,10, 10,15, 15,20), ncol = 2, byrow = TRUE)
#' y_int <- c(0, 2, 5, 8, 10)
#' n_int <- rep(10, 5)
#' result_int <- Griffiths(t_int, y_int, n_int, tau, boot_num = 30)
#' result_int$par
#' result_int$CIs
#'
#' @importFrom stats optim quantile integrate dbinom
#' @export
Griffiths <- function(t, y, n, tau, par_init = c(gamma0 = 0.1, gamma1 = 1), boot_num = 1000) {
  loglik <- function(par, tau, t, y, n) {
    gamma0 <- par[1]
    gamma1 <- par[2]
    if (ncol(t) == 2) { # age buckets
      a <- t[,1]
      b <- t[,2]

      # simplified equation for pi(t) = 1 - e^{inegral foi(t) dt} - see working in goodnotes
      pi_fun <- function(t, gamma0, gamma1, tau) {
        out <- ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau))))
        return(out)
      }

      # we wwant to integrate pi(t) and divide by 1/b-a to get the average pi_t to plug into the binomial distribution (p parameter)
      integrate_pi <- function(a, b, gamma0, gamma1, tau) {
        if (b <= tau) { # then whole age range is below tau, therefore 100% have immunity BUT none of them were previously infected and therefore it makes sense that pi(t) = 0, since pi(t) is the proportion of previously infected (and immune) individuals prior to age t
          return(0)
        } else if (a <= tau && b > tau) { # tau is within the given age bracket
          integrate(pi_fun, lower = tau, upper = b, gamma0 = gamma0, gamma1 = gamma1, tau = tau)$value
        } else {
          integrate(pi_fun, lower = a, upper = b, gamma0 = gamma0, gamma1 = gamma1, tau = tau)$value
        }
      }

      # vector containing average pi_t for each age bucket in t
      pi_t <- mapply(integrate_pi, a, b, MoreArgs = list(gamma0 = gamma0, gamma1 = gamma1, tau = tau)) / (b - a)
    } else { # t is a vector
      pi_t <- pi_t <- 1 - exp((-gamma0*(t^2-tau^2)/2 - gamma0*gamma1*(t-tau))*(t>tau))
    }

    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- dbinom(y, n, pi_t, log=TRUE) # or ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  params <- optim(par=par_init, fn=loglik, tau=tau, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)

  boot_results <- lapply(1:boot_num, function(b) {
    boot_samp <- bootstrap_samples[[b]]

    result <- tryCatch(
      optim(par = par_init, fn = loglik, t = boot_samp$t, tau=tau, y = boot_samp$y, n = boot_samp$n),
      error = function(e) NULL
    )

    if (is.null(result) || result$convergence != 0) {
      return(rep(NA, 2))  # return NA if error occurred or convergence failed
    } else {
      return(result$par)  # return the estimated parameters
    }
  })

  # lapply always outputs a list, while sapply (simplify apply, a wrapped function for lapply that simplifies the list output) simplifies the output into a vector/matrix
  boot_gamma0 <- sapply(boot_results, function(x) x[[1]]) # pull out all the first elements in each list contained in boot_results, and place all the
  boot_gamma1 <- sapply(boot_results, function(x) x[[2]])

  gamma0_CI <- quantile(boot_gamma0, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma1_CI <- quantile(boot_gamma1, probs = c(0.025, 0.975), na.rm = TRUE)

  # to output the fois for the age categories in t, calculate their MLE and confidence intervals
  if (ncol(t) == 1) {
    foi_MLE <- mapply(function(s, gamma0, gamma1, tau) {gamma0*(s+gamma1)*(s>tau)}, s=t, MoreArgs = list(gamma0=params[1], gamma1=params[2], tau=tau))
    #names(foi_MLE) <- as.character(t[,1])
  } else {
    foi_MLE <- mapply(function(s_lower, s_upper, gamma0, gamma1, tau) {1/(s_upper-s_lower) * integrate(function(s) {gamma0 * (s+gamma1)*(s>tau)}, lower = s_lower, upper = s_upper)$value}, s_lower = t[,1], s_upper = t[,2], MoreArgs = list(gamma0 = params[1], gamma1 = params[2], tau = tau))
    #names(foi_MLE) <- apply(t, 1, function(row) paste0("[", row[1], ",", row[2], ")"))
  }

  # foi confidence intervals confidence intervals
  if (ncol(t) == 1) {
    foi_CI <- lapply(t, function(s) {
      boot_foi <- boot_gamma0*(s+boot_gamma1)*as.numeric(s > tau)
      quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)
    }) # note that you don't need to use MoreArgs here since lapply looks at the environment variables, but mapply does not
  } else {
    foi_CI <- lapply(1:nrow(t), function(i) {
      s_lower <- t[i, 1]
      s_upper <- t[i, 2]
      interval_length <- s_upper - s_lower

      # Compute FOI for each bootstrap sample
      boot_foi <- sapply(1:length(boot_gamma0), function(b) {
        integral <- tryCatch(
          integrate(function(s) {boot_gamma0[b] * (s+boot_gamma1[b])*as.numeric(s>tau)}, lower = s_lower, upper = s_upper)$value,
          error = function(e) NA
        )
        return(integral / interval_length)
      })

      quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)
    })
  }


  return(list(par=list(gamma0=params[1], gamma1=params[2], foi=foi_MLE), CIs=list(gamma0_CI=gamma0_CI, gamma1_CI=gamma1_CI, foi_CI=foi_CI), boot_params=list(boot_gamma0=boot_gamma0, boot_gamma1=boot_gamma1)))
}
