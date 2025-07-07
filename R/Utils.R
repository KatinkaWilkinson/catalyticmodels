#' ALTER THIS LATER: Generate Bootstrap Resamples from Count Data
#'
#' Produces a list of bootstrap resamples by sampling with replacement from the pseudo-individual-level data.
#' Each resample is returned in aggregated count form.
#'
#' @param t A numeric vector of ages (or identifiers for age groups).
#' @param y A numeric vector of seropositive counts for each age group.
#' @param n A numeric vector of total sample sizes for each age group.
#' @param num_boot Number of bootstrap resamples to generate.
#'
#' @return A list of length \code{num_boot}. Each element is itself a list with:
#' \describe{
#'   \item{\code{t}}{Vector of unique age values in the bootstrap sample.}
#'   \item{\code{y}}{Bootstrap seropositive counts.}
#'   \item{\code{n}}{Bootstrap total sample sizes.}
#' }
#'
#' @examples
#' t <- c(1, 2, 3)
#' y <- c(2, 3, 4)
#' n <- c(4, 5, 6)
#' create_boot_samps(t, y, n, num_boot = 100)
#'
#' @export
create_boot_samps <- function(t, y, n, num_boot) {
  stopifnot(requireNamespace("future.apply", quietly = TRUE))
  future::plan(future::multisession)

  # Step 1: Build binary data for each age group
  binary_dat <- mapply(function(yi, ni) {
    c(rep(1, yi), rep(0, ni - yi))
  }, y, n, SIMPLIFY = FALSE)

  # Step 2: Bootstrap positive counts for each group
  # One parallel job per group: replicate num_boot times
  group_boots <- future.apply::future_lapply(
    binary_dat,
    function(bin_dat) {
      replicate(num_boot, sum(sample(bin_dat, replace = TRUE)))
    },
    future.seed = TRUE
  )

  # Step 3: Combine bootstraps by sample (i.e., row-wise)
  # Each row becomes one bootstrap dataset
  boot_y_matrix <- do.call(cbind, group_boots)

  # Step 4: Create a list of bootstrap datasets
  boot_list <- lapply(1:num_boot, function(i) {
    list(
      t = t,
      y = boot_y_matrix[i, ],
      n = n
    )
  })

  return(boot_list)
}
create_boot_samps(t, y, n, num_boot)
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
