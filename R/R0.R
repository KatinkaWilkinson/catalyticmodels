R0 <- function(cat_model, age_bins, pop_per_bin) {
  is_spline <- !is.null(cat_model[["spline_pi_t"]])

  # f(a) = proportion of people in age group
  # z(a) = predicted seroprevalence for group (we used mean)
  calc_R0 <- function(z_a) {
    f_a <- pop_per_bin / sum(pop_per_bin)
    1 / (1 - sum(f_a * z_a))
  }

  if (is_spline) { # mean pi (seroprevalence) for splines
    mean_pi <- function(age_bin ,spline_pi_t) {
      l <- age_bin[1]
      u <- age_bin[2]
      integral <- integrate(function(a) min(max(0,predict(spline_pi_t, a)$y),1), lower = l, upper = u)$value
      integral / (l - u)
    }

    seroprev_per_bin <- apply(age_bins, 1, mean_pi, spline_pi_t = cat_model$spline_pi_t)

    R0_Estimate <- calc_R0(seroprev_per_bin)

    bootsamp_mean_pi <- function(y) {
      bootsamp_spline <- smooth.spline(cat_model$t, y / cat_model$n)
      apply(age_bins, 1, mean_pi, bootsamp_spline) # returns vector of mean pi_s for each age group for one boot samp
    }

    boot_seroprev_per_bin <- apply(cat_model$boot_y, 1, bootsamp_mean_pi) # each row is probably one bootsamp's results now

    boot_R0 <- apply(boot_seroprev_per_bin, 2, calc_R0)
    R0_CI <- quantile(boot_R0, probs = c(0.025,0.975), na.rm = TRUE)

    return(list(R0 = R0_Estimate, R0_CI = R0_CI))

  }



  # find mean pi (seroprevalence) for non-splines
  # note that some catalytic models have
  pi_t <- cat_model$pi_t # this is important because some of the built in group_pi functions integrate pi_t
  foi_t <- cat_model$foi_t
  group_pi <- cat_model$group_pi
  mean_pi <- function(age_bin, par) {
    l <- age_bin[1]
    u <- age_bin[2]
    min(max(group_pi(l,u,par),0),1)
  }


  seroprev_per_bin <- apply(age_bins, 1, mean_pi, par = unlist(cat_model$params_MLE))

  R0_MLE <- calc_R0(seroprev_per_bin)

  boot_seroprev_per_bin <- apply(cat_model$bootparams, 1, function(par) apply(age_bins, 1, mean_pi, par=par))

  boot_R0 <- apply(boot_seroprev_per_bin, 2, calc_R0)

  R0_CI <- quantile(boot_R0, probs = c(0.025, 0.975), na.rm = TRUE)

  return(list(R0 = R0_MLE, R0_CI = R0_CI))
}
