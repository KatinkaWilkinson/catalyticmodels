#' Nonparametric FOI Estimation using Isotonic Regression (Keiding's Method)
#'
#' Estimates the force of infection (FOI) nonparametrically using isotonic regression (via the Pool-Adjacent Violators Algorithm),
#' assuming that seroprevalence increases monotonically with age. The FOI is obtained as the derivative of seroprevalence
#' with respect to age, scaled by the susceptible fraction.
#'
#' @param t A numeric vector of exact ages or midpoint values (assumed to be ordered or will be sorted internally).
#' @param y A numeric vector of seropositive counts at each age.
#' @param n A numeric vector of total sample sizes at each age.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{t}}{The age vector corresponding to each FOI estimate (excluding the first point due to differencing).}
#'   \item{\code{foi}}{The pointwise maximum likelihood estimates of the force of infection.}
#'   \item{\code{foi_CI}}{Bootstrap-based 95% confidence intervals for FOI estimates.}
#' }
#'
#' @details
#' This method estimates cumulative seroprevalence \eqn{\pi(t)} using isotonic regression, ensuring a non-decreasing relationship with age.
#' The FOI is computed as:
#' \deqn{
#' \lambda(t) = \frac{d\pi(t)/dt}{1 - \pi(t)}
#' }
#' Bootstrap resampling is used to estimate uncertainty.
#'
#' @export
NonparametricKeidings <- function(t, y, n) {
  ord <- order(t) # pava assumes ts are ordered
  t <- t[ord]
  y <- y[ord]
  n <- n[ord]
  prop_seropos <- y / n
  pi_t <- Iso::pava(prop_seropos) # Estimate cumulative seroprevalence using isotonic regression
  dpi_t_dt <- diff(pi_t)/diff(t)
  MLE_foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)])

  boot_num <- 1000
  boot_foi_t <- matrix(nrow=boot_num, ncol=length(t)-1)
  bootstrap_samples <- create_boot_samps(t, y, n, boot_num)
  for (b in 1:boot_num){
    boot_samp <- bootstrap_samples[[b]]
    bootsamp_t <- boot_samp$t
    bootsamp_y <- boot_samp$y
    bootsamp_n <- boot_samp$n

    bootsamp_params <- tryCatch({
      prop_seropos <- bootsamp_y / bootsamp_n
      pi_t <- pava(prop_seropos) # Estimate cumulative seroprevalence using isotonic regression
      dpi_t_dt <- diff(pi_t)/diff(bootsamp_t)
      foi_t <- dpi_t_dt / (1 - pi_t[-length(pi_t)])
      foi_t
    }, error = function(e) {
      rep(NA, length(t) - 1)
    })

    boot_foi_t[b, ] <- bootsamp_params
  }
  foi_CI <- apply(boot_foi_t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(t = t[-1], foi = MLE_foi_t, foi_CI=foi_CI))
}
