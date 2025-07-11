% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{create_boot_samps}
\alias{create_boot_samps}
\title{Generate Bootstrap Resamples from Grouped Binary Serological Data}
\usage{
create_boot_samps(t, y, n, num_boot)
}
\arguments{
\item{t}{A numeric vector of exact ages, or a matrix with two columns representing lower and upper bounds of age intervals.}

\item{y}{A numeric vector of seropositive counts for each age group or interval.}

\item{n}{A numeric vector of total sample sizes for each age group or interval.}

\item{num_boot}{Integer specifying the number of bootstrap resamples to generate.}
}
\value{
A list of length \code{num_boot}. Each element is itself a list containing:
\describe{
\item{\code{t}}{The original age structure passed to the function — either a vector or a matrix of intervals.}
\item{\code{y}}{Bootstrapped seropositive counts corresponding to each age group or interval.}
\item{\code{n}}{Total sample sizes for each age group or interval (unchanged across bootstraps).}
}
}
\description{
Efficiently generates bootstrap resamples of seropositive counts using grouped binary data derived from count data.
The input \code{t} can either be a vector of exact ages or a matrix representing age intervals (e.g., [a, b] per row).
}
\details{
This function does \strong{not} disaggregate the data (e.g., using Beers, PCLM, or interpolation). All bootstrap samples are created
within the original age groupings — preserving the structure of the input data.

For each age group, a binary vector is constructed representing individual-level serostatus (1 = seropositive, 0 = seronegative),
and bootstrap samples are drawn with replacement \code{num_boot} times. Parallel processing across age groups is done using
\code{future.apply::future_lapply()}.

This method is designed for situations where the age data is already grouped and should remain so throughout inference.
}
\examples{
# Example with exact ages
t_vec <- c(1, 2, 3)
y_vec <- c(2, 3, 4)
n_vec <- c(4, 5, 6)
boot_samples1 <- create_boot_samps(t_vec, y_vec, n_vec, num_boot = 100)

# Example with age intervals
t_mat <- matrix(c(0,5, 5,10, 10,15), ncol = 2, byrow = TRUE)
y <- c(5, 10, 15)
n <- c(10, 20, 25)
boot_samples2 <- create_boot_samps(t_mat, y, n, num_boot = 100)

}
