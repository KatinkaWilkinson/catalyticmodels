#' Synthesize Individual-Level Raw Data from Aggregated Counts
#'
#' Converts grouped count data into pseudo-individual-level data, useful for bootstrapping procedures.
#' For each age value in \code{t}, generates \code{n[i]} rows with corresponding serostatus values.
#'
#' @param t A numeric vector of ages (or age group identifiers).
#' @param y A numeric vector of seropositive counts corresponding to each age group.
#' @param n A numeric vector of total sample sizes corresponding to each age group.
#'
#' @return A matrix with two columns:
#' \describe{
#'   \item{\code{t}}{Age values repeated \code{n[i]} times for each group.}
#'   \item{\code{Sero-result}}{Binary indicator (0 = seronegative, 1 = seropositive).}
#' }
#'
#' @examples
#' t <- c(1, 2, 3)
#' y <- c(1, 2, 3)
#' n <- c(2, 3, 4)
#' synth_data <- synthesise_raw_data(t, y, n)
#'
#' @export
synthesise_raw_data <- function(t, y, n) {
  # synthesise raw data for bootstrapping
  raw_t <- rep(t, n)
  raw_data <- cbind(raw_t, 0)
  colnames(raw_data) <- c("t", "Sero-result") # Sero-result=0=seronegative, Sero-result=1=seropositive
  for (i in 1:length(t)) {
    start_pt <- sum(n[1:(i-1)])+1
    end_pt <- start_pt + y[i] -1
    raw_data[start_pt:end_pt,2] <- rep(1, y[i])
  }
  return(raw_data)
}

#' Aggregate Individual-Level Data into Counts
#'
#' Reconstructs count data (seropositive and total) from pseudo-individual-level data.
#' Useful for re-aggregating bootstrap resamples back into grouped format.
#'
#' @param t A numeric vector of ages or age group identifiers for each individual.
#' @param seroresult A numeric vector of binary serostatus values (0 or 1).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{t}}{Sorted unique age values.}
#'   \item{\code{y}}{Number of seropositives at each age.}
#'   \item{\code{n}}{Total number of individuals at each age.}
#' }
#'
#' @examples
#' t <- c(1,1,2,2,2,3)
#' seroresult <- c(1,0,1,1,0,1)
#' synthesise_count_data(t, seroresult)
#'
#' @export
synthesise_count_data <- function(t, seroresult) {
  agegroups <- sort(unique(t))
  y <- sapply(agegroups, function(tt) sum(seroresult[t==tt]))
  n <- sapply(agegroups, function(tt) sum(t==tt))
  return(list(t=agegroups, y=y, n=n))
}

#' Generate Bootstrap Resamples from Count Data
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
create_boot_samps <- function(t, y, n, num_boot){ # returns the count data from bootstrapping on raw data
  set.seed(123)
  raw_dat <- synthesise_raw_data(t, y, n)
  boot_results <- list()
  for (b in 1:num_boot) {
    ind <- sample(1:nrow(raw_dat), nrow(raw_dat), replace=TRUE)
    boot_raw_dat <- raw_dat[ind,]
    boot_t <- sort(unique(t))
    boot_y <- sapply(boot_t, function(tt) sum(boot_raw_dat[boot_raw_dat[,1] == tt, 2]))
    boot_n <- sapply(boot_t, function(tt) sum(boot_raw_dat[,1] == tt))
    boot_results[[b]] <- list(t=boot_t, y=boot_y, n=boot_n)
  }
  return(boot_results)
}
