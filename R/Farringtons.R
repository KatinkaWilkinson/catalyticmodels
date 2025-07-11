#' Farrington's Time-Dependent Force of Infection Model
#'
#' Implements a non-constant catalytic model where the force of infection (FOI) varies with age following a time-decaying pattern.
#'
#' The FOI is defined as:
#' \deqn{\lambda(t) = (\gamma_0 t + \gamma_1)e^{-\gamma_2 t} + \gamma_1}
#'
#' The cumulative probability of seroconversion up to age \eqn{t} is:
#' \deqn{\pi(t) = 1 - \exp\left(-\int_0^t \lambda(s) ds \right)}
#'
#' For exact ages, the integral is calculated directly. For age intervals \eqn{\[a, b\]}, the average of \eqn{\pi(t)} across the interval is used
#' as the binomial success probability.
#'
#' @param t A matrix of ages. For exact ages, a one-column matrix (e.g., `matrix(age_values, ncol = 1)`). For age intervals, a two-column matrix where each row is a lower and upper bound.
#' @param y A numeric vector of the number of seropositive individuals in each age group or interval.
#' @param n A numeric vector of the number of individuals sampled in each group or interval.
#' @param par_init Named vector of initial parameter values. Defaults to \code{c(gamma0 = 0.1, gamma1 = 1, gamma2 = 0.1)}.
#' @param boot_num Number of bootstrap replicates for computing confidence intervals. Default is 1000.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{par}}{Maximum likelihood estimates of \code{gamma0}, \code{gamma1}, \code{gamma2}, and the FOI for each age group.}
#'   \item{\code{CIs}}{Bootstrap-based 95% confidence intervals for each parameter and the FOI.}
#'   \item{\code{boot_params}}{A list of vectors containing the bootstrap samples for each parameter.}
#' }
#'
#' @details
#' This model is useful when the FOI is expected to decay or level off with age after initially rising.
#' The integration of \eqn{\lambda(t)} over age is done numerically unless a closed-form is available.
#' Confidence intervals are obtained using nonparametric bootstrapping.
#'
#' This function depends on an external utility \code{create_boot_samps()} that generates bootstrap samples in the form:
#' \code{list(list(t = ..., y = ..., n = ...), ...)}.
#'
#' @examples
#' # Exact age example
#' t <- matrix(c(1, 3, 5, 7, 9), ncol = 1)
#' y <- c(0, 2, 5, 8, 10)
#' n <- rep(10, 5)
#' result <- Farringtons(t, y, n, boot_num = 20) # kept low so that example can run quickly
#' result$par
#' result$CIs
#'
#' # Age interval example
#' t_int <- matrix(c(0, 2, 2, 5, 5, 10, 10, 15, 15, 20), ncol = 2, byrow = TRUE)
#' y_int <- c(0, 1, 4, 7, 9)
#' n_int <- rep(10, 5)
#' result_int <- Farringtons(t_int, y_int, n_int, boot_num = 20)
#' result_int$par
#' result_int$CIs
#'
#' @importFrom stats optim quantile integrate
#' @export
Farringtons <- function(t, y, n, par_init = c(gamma0 = 0.1, gamma1 = 1, gamma2 = 0.1), boot_num = 1000) {
  loglik <- function(par, t, y, n) {
    gamma0 <- par[1]
    gamma1 <- par[2]
    gamma2 <- par[3]
    if (ncol(t)==2) {
      a <- t[,1]
      b <- t[,2]

      # function to find integral of the foi from 0 to t - integral solution found online
      pi_fun <- function(t, gamma0, gamma1, gamma2) {
        1-exp(-1*(exp(-gamma2*t)*((gamma1*gamma2^2*t-gamma1*gamma2+gamma0)*exp(gamma2*t)-gamma0*gamma2*t+gamma1*gamma2-gamma0))/(gamma2^2))
      }

      # find the mean of the function over a,b
      pi_t <- mapply(function(a, b) {
        tryCatch(
          integrate(function(t) pi_fun(t, gamma0, gamma1, gamma2), lower = a, upper = b)$value,
          error = function(e) NA
        )
      }, a, b) / (b - a)

    } else {
      pi_t <- sapply(t, function(ti) {
        integrand <- function(s) (gamma0*s - gamma1)*exp(-gamma2*s)+gamma1
        val <- integrate(integrand, 0, ti)$value
        1 - exp(-val)
      })
    }

    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
    return(- sum(ll))
  }

  # MLE
  params <- optim(par=par_init, fn=loglik, t=t, y=y, n=n)$par

  # 95% bootstrap CIs
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)

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
  boot_gamma0 <- sapply(boot_results, function(x) x[[1]]) # pull out all the first elements in each list contained in boot_results, and place all the
  boot_gamma1 <- sapply(boot_results, function(x) x[[2]])
  boot_gamma2 <- sapply(boot_results, function(x) x[[3]])

  gamma0_CI <- quantile(boot_gamma0, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma1_CI <- quantile(boot_gamma1, probs = c(0.025, 0.975), na.rm = TRUE)
  gamma2_CI <- quantile(boot_gamma2, probs = c(0.025, 0.975), na.rm = TRUE)


  # foi MLE for each t

  if (ncol(t) == 1) {
    foi_MLE <- mapply(function(s, gamma0, gamma1, gamma2) {(gamma0*s + gamma1)*exp(-gamma2*s)+gamma1}, s=t, MoreArgs = list(gamma0=params[1], gamma1=params[2], gamma2=params[3]))
    #names(foi_MLE) <- as.character(t[,1])
  } else {
    foi_MLE <- mapply(function(s_lower, s_upper, gamma0, gamma1, gamma2) {1/(s_upper-s_lower) * integrate(function(s) {(gamma0*s + gamma1)*exp(-gamma2*s)+gamma1}, lower = s_lower, upper = s_upper)$value}, s_lower = t[,1], s_upper = t[,2], MoreArgs = list(gamma0 = params[1], gamma1 = params[2], gamma2 = params[3]))
    #names(foi_MLE) <- apply(t, 1, function(row) paste0("[", row[1], ",", row[2], ")"))
  }

  # foi confidence intervals confidence intervals
  if (ncol(t) == 1) {
    foi_CI <- lapply(t, function(s) {
      boot_foi <- (boot_gamma0*s + boot_gamma1)*exp(-boot_gamma2*s)+boot_gamma1
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
          integrate(function(s) {(boot_gamma0[b]*s + boot_gamma1[b])*exp(-boot_gamma2[b]*s)+boot_gamma1[b]}, lower = s_lower, upper = s_upper)$value,
          error = function(e) NA
        )
        return(integral / interval_length)
      })

      # CI from the distribution of bootstrap FOIs
      quantile(boot_foi, probs = c(0.025, 0.975), na.rm = TRUE)
    })
  }

  return(list(par=list(gamma0=params[1], gamma1=params[2], gamma2=params[3], foi=foi_MLE), CIs=list(gamma0_CI=gamma0_CI, gamma1_CI=gamma1_CI, gamma2_CI=gamma2_CI, foi_CI=foi_CI), boot_params=list(boot_gamma0=boot_gamma0, boot_gamma1=boot_gamma1, boot_gamma2=boot_gamma2)))
}
