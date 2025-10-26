#' Generate Bootstrap Resamples from Grouped Binary Serological Data
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
#' @export
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

neg_total_binom_loglik <- function(par, pi_t, group_pi, t, y, n, rho, param_names) {
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

## NEW VERSION BELOW!!! Maybe uncomment!!!

# neg_total_binom_loglik <- function(par, pi_t, group_pi, t, y, n, rho) {
#   # Evaluate pi
#   pi <- tryCatch({
#     if (is.null(dim(t)) || ncol(t) == 1) {
#       sapply(t, function(x) pi_t(x, par))
#     } else {
#       mapply(function(a, b) group_pi(a, b, par), t[,1], t[,2])
#     }
#   }, error = function(e) rep(NA, length(y)))
#
#   # If pi is bad, return penalty
#   if (any(!is.finite(pi)) || any(pi < 0) || any(pi > 1)) return(1e6)
#
#   p_seropos_result <- pmin(pmax(rho * pi, 1e-8), 1 - 1e-8)
#
#   ll <- dbinom(y, size = n, prob = p_seropos_result, log = TRUE)
#
#   total_ll <- sum(ll)
#
#   if (!is.finite(total_ll)) return(1e6)
#
#   return(-total_ll)
# }
