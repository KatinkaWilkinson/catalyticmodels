#' Nonparametric Estimation of Force of Infection via Splines
#'
#' Estimates the force of infection (FOI) as a function of age using a smoothed spline fit
#' to the observed seroprevalence. This approach assumes exact ages and does not yet account
#' for age intervals or grouped data.
#'
#' The force of infection is estimated as:
#' \deqn{
#' \text{FOI}(t) = \frac{d\pi(t)/dt}{1 - \pi(t)}
#' }
#' where \eqn{\pi(t)} is the estimated seroprevalence at age \eqn{t}.
#'
#' @param t A numeric vector of exact ages.
#' @param y A numeric vector of seropositive counts at each age.
#' @param n A numeric vector of total sample sizes at each age.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{t}}{Vector of age points (excluding last grid point).}
#'   \item{\code{foi}}{Estimated force of infection at each age.}
#'   \item{\code{foi_CI}}{Matrix with 95% bootstrap confidence intervals (rows = lower and upper bounds).}
#' }
#'
#' @details
#' The function fits a smoothed spline to \eqn{y/n} and uses numerical differentiation to approximate
#' the derivative of seroprevalence. This is then used to estimate the force of infection.
#'
#' Note: This version assumes input ages represent point values. To extend this to age intervals,
#' one could integrate the spline over each interval or use midpoints as proxies.
#'
#' @export
NonparametricSplines <- function(t, y, n) {
  spline_pi_t <- smooth.spline(t, y/n) # models how seropositive changes with age
  t_grid <- seq(min(t), max(t), length.out = 100) # we only need this if we need the foi over more granular age groups than what we inputted into the function... is this really necessary??? otherwise just use t vector instead of t_grid

  # MLE
  pi_t <- predict(spline_pi_t, t_grid)$y
  dpi_t_dt <- diff(pi_t)/diff(t_grid)
  MLE_foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)])

  # 95% bootstrap CI
  boot_num <- 1000
  boot_foi_t <- matrix(nrow=boot_num, ncol=length(t_grid)-1)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n
    bootsamp_params <- tryCatch({
      spline_pi_t <- smooth.spline(bootsamp_t, bootsamp_y / bootsamp_n)
      pi_t <- predict(spline_pi_t, t_grid)$y
      dpi_t_dt <- diff(pi_t) / diff(t_grid)
      foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)])
      foi_t
    }, error = function(e) {
      rep(NA, length(t_grid) - 1)
    })

    boot_foi_t[b, ] <- bootsamp_params
  }
  foi_CI <- apply(boot_foi_t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(t=t_grid[-length(t_grid)], foi=MLE_foi_t, foi_CI=foi_CI))
}
