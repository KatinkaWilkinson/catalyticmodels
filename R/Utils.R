#' Generate Bootstrap Resamples from Grouped Binary Serological Data (helper)
#'
#' Efficiently generates bootstrap resamples of seropositive counts using grouped binary data derived from count data.
#' The input \code{t} can either be a vector of exact ages or a matrix representing age intervals (e.g., \[a, b\] per row).
#'
#' This function does **not** disaggregate the data (e.g., using Beers, PCLM, or interpolation). All bootstrap samples are created
#' within the original age groupings — preserving the structure of the input data.
#'
#' For each age group, a binary vector is constructed representing individual-level serostatus (1 = seropositive, 0 = seronegative),
#' and bootstrap samples are drawn with replacement \code{num_boot} times. Parallel processing across age groups is done using
#' \code{future.apply::future_lapply()}.
#'
#' @param t A numeric vector of exact ages, or a matrix with two columns representing lower and upper bounds of age intervals.
#' @param y A numeric vector of seropositive counts for each age group or interval.
#' @param n A numeric vector of total sample sizes for each age group or interval.
#' @param num_boot Integer specifying the number of bootstrap resamples to generate.
#'
#' @return A list of length \code{num_boot}. Each element is itself a list containing:
#' \describe{
#'   \item{\code{y}}{Bootstrapped seropositive counts corresponding to each age group or interval.}
#' }
#'
#' @details
#' This method is designed for situations where the age data is already grouped and should remain so throughout inference.
#'
#' @examples
#' # Example with exact ages
#' t_vec <- c(1, 2, 3)
#' y_vec <- c(2, 3, 4)
#' n_vec <- c(4, 5, 6)
#' boot_samples1 <- create_boot_samps(t_vec, y_vec, n_vec, num_boot = 100)
#'
#' # Example with age intervals
#' t_mat <- matrix(c(0,5, 5,10, 10,15), ncol = 2, byrow = TRUE)
#' y <- c(5, 10, 15)
#' n <- c(10, 20, 25)
#' boot_samples2 <- create_boot_samps(t_mat, y, n, num_boot = 100)
#'
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply
#' @noRd
create_boot_samps <- function(t, y, n, num_boot) {
  stopifnot(requireNamespace("future.apply", quietly = TRUE))
  future::plan(future::multisession)

  # binary_dat is a list where each entry is a binary vector of data for one age group
  binary_dat <- mapply(function(yi, ni) {c(rep(1, yi), rep(0, ni - yi))}, y, n, SIMPLIFY = FALSE)

  # For each element in the binary_dat list (representing the binary sero data for a single age group/category), we draw num_boot bootstrap samples with replacement and sum the values in each sample to obtain num_boot bootstrapped seropositive counts. future_lapply creates a different thread for each element in the list (each age category), so that each thread is working on generating num_boot bootstrapped y values for one of the age categories. Group_boots is now a list where each element is a vector of num_boot y values
  group_boots <- future.apply::future_lapply(binary_dat,
    function(bin_dat) {replicate(num_boot, sum(sample(bin_dat, replace = TRUE)))},
    future.seed = TRUE # for safe parallel use of random numbers in sample() func
  )

  # Combine the list of bootstrap vectors into a matrix: rows correspond to bootstrap iterations, columns to age groups
  boot_y_matrix <- do.call(cbind, group_boots)

  # We want the output of the function to be a list where each element in the list is itself a list containing the t, y, n triplet for one bootstrap sample.
  boot_list <- lapply(1:num_boot, function(i) {list(y = boot_y_matrix[i, ])})

  future::plan(future::sequential) # close plan

  return(boot_list)
}

#' Negative total binomial log-likelihood (helper)
#'
#' Computes the negative log-likelihood for binomial seroprevalence counts,
#' given a prevalence model \code{pi_t()} (for exact-age data) or
#' \code{group_pi()} (for age-interval data). Used internally by
#' \code{FoiFromCatalyticModel()} during maximum-likelihood estimation.
#'
#' @param par Named numeric vector of model parameters. If names are missing,
#'   they are assigned from \code{param_names}.
#' @param pi_t Function \code{pi_t(t, par)} returning point prevalence at age \code{t}.
#' @param group_pi Function \code{group_pi(a, b, par)} returning average prevalence on \eqn{[a, b)}.
#' @param group_foi Unused here (kept for a consistent signature with other helpers).
#' @param t Numeric vector of ages, or a two-column matrix of age intervals \eqn{[a, b)}.
#' @param y Integer vector of observed seropositive counts per age (or age group).
#' @param n Integer vector of totals tested per age (or age group).
#' @param rho Numeric scalar or \code{NA}. If \code{NA}, an estimated \code{rho}
#'   is taken from \code{par["rho"]}; otherwise the supplied scalar is used.
#' @param param_names Character vector of parameter names to use if \code{par}
#'   arrives unnamed (e.g., from \code{optim}).
#'
#' @details
#' The routine:
#' \enumerate{
#'   \item Builds \code{pi} via \code{pi_t()} for exact-age input, or
#'         via \code{group_pi()} for interval input.
#'   \item Applies the test-adjustment multiplicatively: \eqn{p = \rho \cdot \pi}.
#'   \item Bounds \eqn{p} to \eqn{[1e-8, 1-1e-8]} for numerical stability.
#'   \item Returns \code{-sum(dbinom(y, n, p, log = TRUE))}.
#' }
#' If the total log-likelihood is non-finite, a large penalty (\code{1e6}) is returned.
#'
#' @return Numeric scalar: the negative log-likelihood.
#'
#' @keywords internal
#' @importFrom stats dbinom
#' @noRd
neg_total_binom_loglik <- function(par, pi_t, group_pi, group_foi, t, y, n, rho, param_names) {
  if (is.null(names(par))) {names(par) <- param_names} # putting this in because optim drops the names!

  # Compute pi depending on whether t is intervals or points
  if (is.null(dim(t)) || ncol(t) == 1) {
    pi <- sapply(t, function(x) pi_t(x, par))
  } else {
    pi <- mapply(function(a, b) group_pi(a, b, par), t[,1], t[,2])
  }

  # Determine rho estimate
  if (is.na(rho)) {
    rho_estimate <- par["rho"]
  } else {
    rho_estimate <- rho
  }

  # Calculate seropositivity probability with rho adjustment
  p_seropos_result <- rho_estimate * pi

  # Bound probabilities away from 0 and 1 for numerical stability
  p_seropos_result <- pmin(pmax(p_seropos_result, 1e-8), 1 - 1e-8)

  # Log-likelihood of observed seropositive counts under binomial

  # LOGLIK VERSION
  ll <- dbinom(y, size = n, prob = p_seropos_result, log = TRUE)

  total_ll <- sum(ll)

  # Penalize non-finite likelihoods to avoid optimization errors
  if (!is.finite(total_ll)) {
    return(1e6)
  }

  # Return negative log-likelihood for minimization
  return(-total_ll)

  # # DEVIANCE VERSION
  # loglik_deviance <- -2*(sum(y*log(p_seropos_result) + (n-y)*log(1-p_seropos_result) - y*log(y/n) - (n-y)*log(1-y/n)))
  #
  # if (!is.finite(loglik_deviance)) {
  #   return(1e6)
  # }
  #
  # return(loglik_deviance)
}


