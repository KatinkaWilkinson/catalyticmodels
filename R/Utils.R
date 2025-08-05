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
#'   \item{\code{t}}{The original age structure passed to the function — either a vector or a matrix of intervals.}
#'   \item{\code{y}}{Bootstrapped seropositive counts corresponding to each age group or interval.}
#'   \item{\code{n}}{Total sample sizes for each age group or interval (unchanged across bootstraps).}
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

neg_total_binom_loglik <- function(par, pi_t, group_pi, t, y, n, rho) {
  if (is.null(dim(t)) || ncol(t) == 1) {
    pi <- sapply(t, function(x) pi_t(x, par))
  } else {
    pi <- mapply(function(a, b) group_pi(a, b, par), t[,1], t[,2])
  }

  p_seropos_result <- rho * pi
  p_seropos_result <- pmin(pmax(p_seropos_result, 1e-8), 1 - 1e-8)

  ll <- dbinom(y, size = n, prob = p_seropos_result, log = TRUE)

  total_ll <- sum(ll)

  if (!is.finite(total_ll)) {
    return(1e6)  # Penalize: bad parameter set
  }

  return(-total_ll)
}


#
# create_boot_samps_old_and_incorrect <- function(t, y, n, num_boot){
#   set.seed(123)
#   raw_dat <- synthesise_raw_data(t, y, n)
#   boot_results <- list()
#
#   for (b in 1:num_boot) {
#     ind <- sample(1:nrow(raw_dat), nrow(raw_dat), replace = TRUE)
#     boot_raw_dat <- raw_dat[ind, ]
#
#     if (is.matrix(t)) {
#       # Recreate the same intervals
#       age_intervals <- t
#       boot_y <- numeric(nrow(age_intervals))
#       boot_n <- numeric(nrow(age_intervals))
#       for (i in 1:nrow(age_intervals)) {
#         lower <- age_intervals[i, 1]
#         upper <- age_intervals[i, 2]
#         in_interval <- boot_raw_dat$t >= lower & boot_raw_dat$t < upper
#         boot_y[i] <- sum(boot_raw_dat[in_interval,2])
#         boot_n[i] <- sum(in_interval)
#       }
#       boot_results[[b]] <- list(t = age_intervals, y = boot_y, n = boot_n)
#
#     } else {
#       # For exact ages
#       boot_t <- t
#       boot_y <- sapply(boot_t, function(tt) sum(boot_raw_dat[boot_raw_dat[,1] == tt, 2]))
#       boot_n <- sapply(boot_t, function(tt) sum(boot_raw_dat[,1] == tt))
#       boot_results[[b]] <- list(t = boot_t, y = boot_y, n = boot_n)
#     }
#   }
#
#   return(boot_results)
# }

# synthesise_raw_data <- function(t, y, n) {
#   total_people <- sum(n)
#   raw_mat <- matrix(0, nrow = total_people, ncol = 2)  # columns: age, seroresult
#
#   if (is.matrix(t) && ncol(t) == 2) {
#     # Interval data
#     current_index <- 1
#     for (i in 1:nrow(t)) {
#       num_people <- n[i]
#       num_positive <- y[i]
#       num_negative <- num_people - num_positive
#
#       ages <- runif(num_people, min = t[i, 1], max = t[i, 2])
#       sero_results <- sample(c(rep(1, num_positive), rep(0, num_negative)))
#
#       rows <- current_index:(current_index + num_people - 1)
#       raw_mat[rows, 1] <- ages
#       raw_mat[rows, 2] <- sero_results
#
#       current_index <- current_index + num_people
#     }
#
#   } else {
#     # Exact ages
#     raw_t <- rep(t, n)
#     sero_results <- numeric(length(raw_t))  # preallocate
#
#     index <- 1
#     for (i in seq_along(t)) {
#       num_people <- n[i]
#       num_positive <- y[i]
#       num_negative <- num_people - num_positive
#
#       sero_vals <- sample(c(rep(1, num_positive), rep(0, num_negative)))
#
#       sero_results[index:(index + num_people - 1)] <- sero_vals
#       index <- index + num_people
#     }
#
#     raw_mat[, 1] <- raw_t
#     raw_mat[, 2] <- sero_results
#   }
#
#   # Final checks before returning
#   raw_data <- data.frame(t = raw_mat[, 1], seroresult = raw_mat[, 2])
#
#   # Debug print
#   print(table(raw_data$t))
#   print(table(raw_data$seroresult))
#   return(raw_data)
# }

# synthesise_count_data <- function(t, seroresult) {
#   agegroups <- sort(unique(t))
#   y <- sapply(agegroups, function(tt) sum(seroresult[t==tt]))
#   n <- sapply(agegroups, function(tt) sum(t==tt))
#   return(list(t=agegroups, y=y, n=n))
# }
