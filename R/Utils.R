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

synthesise_count_data <- function(t, seroresult) {
  agegroups <- sort(unique(t))
  y <- sapply(agegroups, function(tt) sum(seroresult[t==tt]))
  n <- sapply(agegroups, function(tt) sum(t==tt))
  return(list(t=agegroups, y=y, n=n))
}

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