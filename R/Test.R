# ------------------------------------------------------------------------------------
# test exact ages with all integrals provided
# this example uses Griffith's

t <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n <- rep(50, length(t))
y <- c(0, 0, 1, 2, 5, 10, 20, 30, 40, 45)

tau <- 2  # maternal antibody cutoff

# Griffiths pi_t function
pi_t <- function(t, par) {
  tau <- 2 # in this case!
  gamma0 <- par[1]
  gamma1 <- par[2]
  return(ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau)))))
}

foi_t <- function(t, par) {
  tau <- 2 # in this case!
  gamma0 <- par[1]
  gamma1 <- par[2]
  return(ifelse(t <= tau, 0, gamma0 * (t + gamma1)))
}

par_init <- c(gamma0=0.1, gamma1=0.5) # note that whatever you call your pars in par_init will be what the params_MLE are labelled as (I did not code it to be this way, but it happens yay)

result <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  pi_t = pi_t,
  foi_t = foi_t,
  par_init = par_init
)

result$foi_CI

# ------------------------------------------------------------------------------------
# test age buckets with all integrals provided
# this example uses Griffith's

t <- matrix(c(1,4, 4,8, 8,11, 11,15, 15,18, 18,22, 22,25, 25,29, 29,32, 32,35), ncol=2, byrow = TRUE)
n <- rep(50, length(t))
y <- c(0, 0, 1, 2, 5, 10, 20, 30, 40, 45)

# method of obtaining this function: Put griffith's equation into the integral calculator (google "integral calculator"). copied the latex of the definite integral from a to b and put it into chatGPT. Asked chat to turn it into a function.

group_pi <- function(a, b, par) {
  gamma0 <- par[1]
  gamma1 <- par[2]
  tau <- 2

  if (b <= tau) return(0)

  # Adjust bounds if partially below tau
  a_adj <- max(a, tau)
  b_adj <- b

  sqrt_pi <- sqrt(pi)
  sqrt2 <- sqrt(2)
  sqrt_gamma0 <- sqrt(gamma0)

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  exp_term <- exp((gamma0 * a_adj^2) / 2 + gamma0 * gamma1 * a_adj + (gamma0 * gamma1^2) / 2)

  erf1 <- erf(sqrt_gamma0 * (sqrt2 * a_adj + sqrt2 * gamma1) / 2)
  erf2 <- erf(sqrt_gamma0 * (sqrt2 * b_adj + sqrt2 * gamma1) / 2)

  A_term <- sqrt_pi * sqrt_gamma0 * sqrt2 * exp_term * (erf1 - erf2)
  B_term <- 2 * gamma0 * (b_adj - a_adj)

  integral <- (A_term + B_term) / (2 * gamma0)

  avg_pi <- integral / (b - a)

  return(avg_pi)
}

group_foi <- function(a, b, par) {
  gamma0 <- par[1]
  gamma1 <- par[2]
  tau <- 2

  if (b <= tau) return(0)

  # Adjust bounds if partially below tau
  a_adj <- max(a, tau)
  b_adj <- b
  integral <- ((b_adj - a_adj) * (2 * gamma0 * gamma1 + (b_adj + a_adj) * gamma0)) / 2
  avg_foi <- integral / (b - a)

  return(avg_foi)
}

par_init <- c(gamma0=0.1, gamma1=0.5) # note that whatever you call your pars in par_init will be what the params_MLE are labelled as (I did not code it to be this way, but it happens yay)

result <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  group_pi = group_pi,
  group_foi = group_foi,
  par_init = par_init
)

result

# ----------------------------------------------------------------------------------------

# TEST type = something
# exact ages
t <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
n <- rep(50, length(t))
y <- c(3, 5, 15, 20, 25, 35, 35, 39, 40, 45)

result_muenchgeneral <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "MuenchGeneral",
  model_fixed_params = list()
)

result_muenchgeneral

plot_foi_grid(result_muenchgeneral, 1, 10, confint=TRUE)

result_muenchresticted <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "MuenchRestricted",
  lower = 0.001,
  upper = 3,
  model_fixed_params = list(k=1, l=1)
)

result_muenchresticted

plot_foi_grid(result_muenchresticted, 1, 10, confint=TRUE)


result_griffiths <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "Griffiths",
  model_fixed_params = list(tau=2)
)

result_griffiths

plot_foi_grid(result_griffiths, 1, 10, confint=TRUE)

result_farringtons <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "Farringtons"
)

result_farringtons

plot_foi_grid(result_farringtons, 1, 10, confint=TRUE)

result_PiecewiseConstant <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "PiecewiseConstant",
  model_fixed_params = list(upper_cutoffs = c(3, 7, Inf))
)

result_PiecewiseConstant

plot_foi_grid(result_PiecewiseConstant, 1, 10, confint=TRUE)

result_splines <- FoiFromCatalyticModel(
  t=t,
  y=y,
  n=n,
  type = "Splines"
)

result_splines

plot_foi_grid(result_splines, 1, 10, confint=TRUE)


plot(result_splines$foi_grid ~ seq(min(t), max(t), length.out = 100),
     type = "l", col = "blue", lwd = 2,
     ylab = "Force of Infection", xlab = "Age",
     main = "FOI from Splines (Exact Age Data)")

plot_foi_grid(list(Muench = result_muenchgeneral,Griffiths = result_griffiths,Farringtons=result_farringtons,`Piecewise Constant` = result_PiecewiseConstant), from = 1, to = 10, confint = TRUE)

plot_foi_grid(list(`Muench General` = result_muenchgeneral, `Muench Restricted` = result_muenchresticted, Griffiths = result_griffiths,`Piecewise Constant` = result_PiecewiseConstant), from = 1, to = 10, confint = FALSE)


# age bands
t <- matrix(c(1,4, 4,8, 8,11, 11,15, 15,18, 18,22, 22,25, 25,29, 29,32, 32,35), ncol=2, byrow = TRUE)
n <- rep(50, nrow(t))
y <- c(0, 0, 1, 2, 5, 10, 20, 30, 40, 45)

result_muenchgeneral_agebands <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "MuenchGeneral",
  model_fixed_params = list()
)

result_muenchgeneral_agebands

result_muenchresticted_agebands <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "MuenchRestricted",
  lower = 0.001,
  upper = 3,
  model_fixed_params = list(k=1, l=1)
)

result_muenchresticted_agebands

result_griffiths_agebands <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "Griffiths",
  model_fixed_params = list(tau=2),
  lower = c(1e-6, -20),
  upper = c(0.5, 10)
)

result_griffiths_agebands

result_farringtons_agebands <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "Farringtons",
  lower = c(gamma0 = 1e-6, gamma1 = 1e-6, gamma2 = 1e-3),
  upper = c(gamma0 = 10, gamma1 = 10, gamma2 = 2),
  boot_num = 200 # because of the integrate function, this one takes looong!
)

result_farringtons_agebands # look into this!! What should gamma0,gamma1,gamma2 restrictions be?

result_PiecewiseConstant_agebands <- FoiFromCatalyticModel(
  t = t,
  y = y,
  n = n,
  type = "PiecewiseConstant",
  model_fixed_params = list(upper_cutoffs = c(13, 25, Inf)),
  lower = c(-1,-1,-1),
  upper = c(5,5,5),
  maxit = 1000,
  factr = 1e6
)

result_PiecewiseConstant_agebands

result_splines_agebands <- FoiFromCatalyticModel(
  t=t,
  y=y,
  n=n,
  type = "Splines"
)

result_splines_agebands

plot(result_splines_agebands$foi_grid ~ seq(min(t), max(t), length.out = 100),
     type = "l", col = "blue", lwd = 2,
     ylab = "Force of Infection", xlab = "Age",
     main = "FOI from Splines (Exact Age Data)")

plot_foi_grid(list(Muench = result_muenchgeneral_agebands,Griffiths = result_griffiths_agebands,Farringtons=result_farringtons_agebands,`Piecewise Constant` = result_PiecewiseConstant_agebands), from = 1, to = 35, confint = TRUE)

plot_foi_grid(list(`Muench General` = result_muenchgeneral_agebands, `Muench Restricted` = result_muenchresticted_agebands, Griffiths = result_griffiths_agebands,`Piecewise Constant` = result_PiecewiseConstant_agebands), from = 1, to = 35, confint = FALSE)


#### SPEED TESTS


# Load required package
library(microbenchmark)

# Closed-form integral using pnorm (i.e., Griffiths' model)
griffiths_integral <- function(a, b, gamma0, gamma1) {
  sqrt_pi <- sqrt(pi)
  sqrt2 <- sqrt(2)
  sqrt_gamma0 <- sqrt(gamma0)

  # Error function via pnorm
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  term1 <- erf(sqrt_gamma0 * (sqrt2 * a + sqrt2 * gamma1) / 2)
  term2 <- erf(sqrt_gamma0 * (sqrt2 * gamma1 + sqrt2 * b) / 2)

  A <- sqrt_pi * sqrt_gamma0 * (
    sqrt2 * exp((gamma0 * a^2)/2 + gamma0 * gamma1 * a + (gamma0 * gamma1^2)/2) * term1 -
      sqrt2 * exp((gamma0 * a^2)/2 + gamma0 * gamma1 * a + (gamma0 * gamma1^2)/2) * term2
  )

  B <- -2 * gamma0 * a + 2 * gamma0 * b

  result <- (A + B) / (2 * gamma0)
  return(result)
}

# Numerical version using base R integrate()
griffiths_fun <- function(t, gamma0, gamma1, tau) {
  ifelse(t <= tau, 0, 1 - exp(-((gamma0 / 2) * (t^2 - tau^2) + gamma0 * gamma1 * (t - tau))))
}

numerical_integral <- function(b, gamma0, gamma1, tau) {
  if (b <= tau) return(0)
  integrate(griffiths_fun, lower = tau, upper = b, gamma0 = gamma0, gamma1 = gamma1, tau = tau)$value
}

# Example parameters
gamma0 <- 0.3
gamma1 <- 2
tau <- 1
b <- 5

# Check output values are close
cat("Closed-form:", griffiths_integral(tau, b, gamma0, gamma1), "\n")
cat("Numerical:  ", numerical_integral(b, gamma0, gamma1, tau), "\n")

# Time them both
microbenchmark(
  closed_form = griffiths_integral(b, gamma0, gamma1, tau),
  numerical = numerical_integral(b, gamma0, gamma1, tau),
  times = 100L
)
