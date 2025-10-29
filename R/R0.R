#' Estimate R0 from a fitted catalytic model by age-bin averaging
#'
#' Computes a crude basic reproduction number \eqn{R_0} from a fitted catalytic
#' model by averaging the model-predicted seroprevalence within user-defined
#' age bins and combining these with population weights. The estimator used is
#' \deqn{R_0 = 1 / \{ 1 - \sum_a f_a \, \bar{\pi}(a) \},}
#' where \eqn{f_a} is the fraction of the population in age bin \eqn{a}, and
#' \eqn{\bar{\pi}(a)} is the mean predicted seroprevalence in that bin.
#'
#' The function supports two classes of fitted objects:
#' \itemize{
#'   \item \strong{Spline models}: identified by \code{cat_model$spline_pi_t}.
#'     Mean seroprevalence per bin is computed by numerical integration of
#'     the spline prediction. Bootstrap uncertainty is propagated by refitting
#'     a spline to each bootstrap sample in \code{cat_model$boot_y} and repeating
#'     the bin means and \(R_0\) calculation.
#'   \item \strong{Closed-form models}: identified by the presence of
#'     \code{cat_model$group_pi}, \code{cat_model$params_MLE}, and
#'     \code{cat_model$bootparams}. Mean seroprevalence per bin is obtained via
#'     \code{group_pi(lower, upper, par)}, then \(R_0\) is computed for the MLE
#'     and for each bootstrap parameter set.
#' }
#'
#' @param cat_model A fitted catalytic model object from your workflow. For
#'   spline models it must provide at least:
#'   \code{$spline_pi_t} (a \code{smooth.spline} for \(\pi(t)\)),
#'   \code{$t} (observation times; vector or two-column matrix of lower/upper
#'   bounds), \code{$n} (trial sizes), and \code{$boot_y} (bootstrap counts laid
#'   out row-wise for refitting). For closed-form models it must provide:
#'   \code{$group_pi} (function \code{group_pi(a,b,par)} returning mean
#'   seroprevalence in \([a,b]\)), \code{$params_MLE} (named numeric),
#'   and \code{$bootparams} (matrix of bootstrap parameter draws; one set per row).
#' @param age_bins Numeric matrix with two columns \code{[lower, upper]} giving
#'   the age-bin edges (same units as the model’s \code{t}). One row per age bin.
#' @param pop_per_bin Numeric vector of population counts per age bin. Length
#'   must equal \code{nrow(age_bins)}. Values are internally normalised to
#'   fractions \eqn{f_a}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{R0}}{Point estimate of \(R_0\).}
#'   \item{\code{R0_CI}}{A length-2 named numeric vector with the 2.5% and 97.5%
#'   bootstrap quantiles.}
#' }
#'
#' @details
#' \itemize{
#'   \item For spline models, bin means \(\bar{\pi}(a)\) are computed as
#'   \(\frac{1}{u-l}\int_l^u \hat{\pi}(t)\,dt\) using \code{stats::integrate} on
#'   the spline’s \code{predict()}ed values. For each bootstrap \code{y} the code
#'   refits a \code{smooth.spline} to either the midpoints (interval data) or
#'   the raw times \code{t} (point data).
#'   \item For closed-form models, \code{group_pi()} is used per bin and the
#'   result is clipped to \code{[0,1]} for numerical safety.
#'   \item The \(R_0\) identity assumes endemic equilibrium and life-long immunity;
#'   interpret with care when these assumptions do not hold.
#' }
#'
#' @section Requirements and assumptions:
#' \itemize{
#'   \item \code{age_bins} has two columns \code{[lower, upper]} with
#'   \code{lower < upper} for all rows.
#'   \item \code{pop_per_bin} is non-negative with at least one positive entry.
#'   \item For spline models, \code{spline_pi_t} is a \code{smooth.spline} whose
#'   \code{predict()} returns a list with element \code{$y} on numeric input.
#' }
#'
#' @examples
#' \dontrun{
#' # Define bins and population
#' age_bins <- matrix(c(0,5, 5,10, 10,15, 15,20), ncol = 2, byrow = TRUE)
#' pop      <- c(1200, 1300, 1100, 1000)
#'
#' # Closed-form example
#' R0(cm_parametric, age_bins, pop)
#'
#' # Spline example
#' R0(cm_spline, age_bins, pop)
#' }
#'
#' @seealso \code{\link[stats]{integrate}}, \code{\link[stats]{smooth.spline}},
#'   \code{\link[stats]{predict}}, \code{\link[stats]{quantile}}
#'
#' @importFrom stats integrate smooth.spline predict quantile
#' @export
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
      pi_t <- function(t, spline_pi_t) {
        predict(spline_pi_t, t)$y
      }

      group_pi <- function(a, b, spline_pi_t) {
        integral <- integrate(pi_t, a, b, spline_pi_t=spline_pi_t)$value
        integral/(b-a)
      }
      group_pi(l, u, spline_pi_t)
    }

    seroprev_per_bin <- apply(age_bins, 1, mean_pi, spline_pi_t = cat_model$spline_pi_t)

    R0_Estimate <- calc_R0(seroprev_per_bin)

    bootsamp_mean_pi <- function(y) {
      if (is.matrix(cat_model$t) && ncol(cat_model$t) == 2) {
        bootsamp_spline <- smooth.spline((cat_model$t[,1] + cat_model$t[,2])/2 , y / cat_model$n)
      } else (
        bootsamp_spline <- smooth.spline(cat_model$t , y / cat_model$n)
      )
      apply(age_bins, 1, mean_pi, bootsamp_spline) # returns vector of mean pi_s for each age group for one boot samp
    }

    boot_y <- matrix(cat_model$boot_y, ncol = length(cat_model$n), byrow = TRUE)
    boot_seroprev_per_bin <- apply(boot_y, 1, bootsamp_mean_pi) # each row is probably one bootsamp's results now

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
