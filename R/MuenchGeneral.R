#' Muench's Catalytic Model
#'
#' Fits the catalytic model: π(t) = k(l - e^{-λt})
#'
#' @param t Vector of age groups
#' @param y Seropositive counts
#' @param n Sample sizes
#' @return A list with MLE estimates, CIs, and bootstrap samples
#' @export
MuenchGeneral <- function(t, y, n) {
  # t = (t1, ..., tN) = vector of N age groups
  # y = (y1, ..., yN) = vector of number of seropositive in each age group
  # n = (n1, ..., nN) = vector of total counts per age group
  loglik <- function(par, t, y, n) {
    k <- par[1]
    l <- par[2]
    foi <- par[3]
    pi_t <- k * (l - exp(-foi*t))
    pi_t <- pmin(pmax(pi_t, 1e-8), 1 - 1e-8) # to ensure that you don't get a pi_t of 0 or 1 - this avoids a log(0) error
    ll <- y*log(pi_t) + (n-y)*log(1-pi_t)
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
