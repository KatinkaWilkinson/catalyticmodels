#' Vynnycky Catalytic Model B (Piecewise Constant FOI, Assay Sensitivity = 1)
#'
#' Fits a piecewise constant force of infection (FOI) catalytic model assuming full assay sensitivity (\eqn{ρ = 1}).
#'
#' The seroprevalence function is modeled as:
#' \deqn{
#' \pi(t) =
#' \begin{cases}
#'   1 - \exp(-\lambda_y \cdot (t - 0.5)), & \text{for } t < 13 \\
#'   1 - \exp(-\lambda_y \cdot 12.5 - \lambda_o \cdot (t - 13)), & \text{for } t \geq 13
#' \end{cases}
#' }
#'
#' This structure captures different FOI values for children (\eqn{\lambda_y}) and older individuals (\eqn{\lambda_o}).
#'
#' @param t A numeric vector of age groups.
#' @param y A numeric vector of seropositive counts corresponding to each age group.
#' @param n A numeric vector of total sample sizes corresponding to each age group.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{par}}{Named list of maximum likelihood estimates: \code{foi_y}, \code{foi_o}.}
#'   \item{\code{CIs}}{Named list of 95% bootstrap confidence intervals for each parameter.}
#'   \item{\code{boot_params}}{Named list containing the bootstrap samples for \code{foi_y} and \code{foi_o}.}
#' }
#'
#' @export
VynnyckyB <- function(t, y, n) {
  loglik <- function(par, t, y, n) {
    foi_y <- par[1]
    foi_o <- par[2]
    pi_t <- 1-exp(-foi_y*pmin(t-0.5, 12.5) -foi_o*pmax(0,t-13))
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
