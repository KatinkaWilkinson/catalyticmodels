% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MuenchRestricted.R
\name{MuenchRestricted}
\alias{MuenchRestricted}
\title{Muench's Catalytic Model with Constant Force of Infection (Restricted)}
\usage{
MuenchRestricted(
  t,
  y,
  n,
  par_init = c(foi = 0.1),
  search_range = c(lower = -1, upper = 3)
)
}
\arguments{
\item{t}{A numeric vector of exact ages, or a matrix with two columns specifying the lower and upper bounds of age intervals.}

\item{y}{A numeric vector of seropositive counts corresponding to each age or age interval.}

\item{n}{A numeric vector of total sample sizes for each age or age interval.}

\item{par_init}{Optional initial value for \code{foi} (default is 0.1).}

\item{search_range}{A numeric vector of length 2 giving the lower and upper bounds for optimization via Brent's method (default is c(-1, 3)).}
}
\value{
A list with:
\describe{
\item{\code{par}}{Maximum likelihood estimate of the force of infection (\code{foi}).}
\item{\code{CIs}}{Bootstrap-based 95\% confidence interval for \code{foi}.}
\item{\code{boot_params}}{Vector of bootstrap estimates for \code{foi}.}
}
}
\description{
Fits the classical Muench catalytic model under the assumption of a constant force of infection (FOI),
where \eqn{\pi(t) = 1 - \exp(-\lambda t)}. This formulation assumes full susceptibility from birth
and that infection accumulates exponentially with age until saturation.
}
\details{
For age intervals \eqn{[a, b]}, the model estimates the average seroprevalence within the interval as:
\deqn{\bar{\pi}_{[a,b]} = 1 - \frac{e^{-\lambda a} - e^{-\lambda b}}{\lambda(b - a)}}
This is derived analytically to avoid computationally expensive numerical integration.

This restricted form of the Muench model sets the scaling constants \eqn{k = 1} and \eqn{l = 1}, meaning
all individuals are assumed to be fully susceptible and eventually infected if exposed long enough.
This makes it appropriate for endemic infections with no significant maternal immunity or heterogeneity in susceptibility.
Brent’s method is used for 1D optimization for numerical stability and performance.
}
\examples{
# Using exact ages
t <- matrix(1:10, ncol = 1)
y <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
n <- rep(10, 10)
result <- MuenchRestricted(t, y, n)
result$par

# Using age intervals
t_int <- matrix(c(0,1, 1,2, 2,3, 3,4, 4,5), ncol = 2, byrow = TRUE)
y <- c(0, 1, 2, 3, 4)
n <- rep(10, 5)
result2 <- MuenchRestricted(t_int, y, n)
result2$CIs

}
