# Clear the workspace
rm(list=ls())

###############################################
### Example script to use the run_mcmc function
###############################################

# Load necessary libraries
library(coda)       # For trace plots and density plots
library(truncnorm)  # For truncated normal distribution
library(MASS)       # For multivariate normal distribution
library(Matrix)     # For Kronecker product
library(mvtnorm)    # For multivariate normal distribution
library(R.utils)    # For saving files

# Define your data and parameters
n <- 44  # Number of observations
n.iter <- 1000  # Number of MCMC iterations
file_name <- "MCMC_Output"  # Base name for saving files
save_dir <- "/Users/pwill/Dropbox/GitHub/Wood-Duck-Mercury-Contamination/R_Output"  # Directory to save output

# Simulate some data
set.seed(123)  # For reproducibility

# Matrix building materials
m.tmp1 <- matrix(0, 6, 6)
m.tmp1[c(1, 8, 15)] <- 1
m.tmp2 <- matrix(0, 6, 6)
m.tmp2[c(22, 29, 36)] <- 1
m.tmp3 <- matrix(0, 6, 6)
m.tmp3[c(2, 3, 7, 9, 13, 14)] <- 1
m.tmp4 <- matrix(0, 6, 6)
m.tmp4[c(23, 24, 28, 30, 34, 35)] <- 1
m.tmp5 <- matrix(0, 6, 6)
m.tmp5[c(4, 5, 6, 10, 11, 12, 16, 17, 18, 19, 20, 21, 25, 26, 27, 31, 32, 33)] <- 1

# True parameter values for simulation
s.l.truth <- 800
s.b.truth <- sqrt(2.5)
rho.truth <- c(0.7, 0.5, 0.3)
block <- s.l.truth^2 * m.tmp1 +
  s.b.truth^2 * m.tmp2 +
  rho.truth[1] * m.tmp3 * s.l.truth^2 +
  rho.truth[2] * m.tmp4 * s.b.truth^2 +
  rho.truth[3] * s.l.truth * s.b.truth * m.tmp5
Sigma <- kronecker(diag(n), block)

# Generate synthetic data
n.mall <- 5
liver <- rep(rep(1:0, each = 3), n)
species <- c(rep(0, (n - n.mall) * 6), rep(1, n.mall * 6))
flank1 <- rnorm(n, 8517.324, 8709.865)
flank2 <- rnorm(n, 10013.11, 11485.8)
axillary1 <- c(rnorm(n - n.mall, 13958.39, 16311.44), rep(0, n.mall))
axillary2 <- c(rnorm(n - n.mall, 13213.99, 15451.69), rep(0, n.mall))
breast1 <- rnorm(n, 9873, 12384)
breast2 <- rnorm(n, 9873, 12207)
flank <- rep((flank1 + flank2) / 2, each = 6)
axillary <- rep((axillary1 + axillary2) / 2, each = 6)
breast <- rep((breast1 + breast2) / 2, each = 6)
X <- matrix(1, length(breast), 1)
X <- cbind(X,
           liver,
           scale(flank),
           scale(breast),
           species,
           species * liver,
           species * scale(flank),
           (1 - species) * scale(axillary),
           species * scale(breast))
beta.truth <- rnorm(ncol(X), 0, 2)
mu <- matrix(X %*% beta.truth, n * 6, 1)
y <- matrix(rmvnorm(1, mean = mu, Sigma), , 1)

# Set initial values, tuning parameters, and prior parameters
initial_values <- list(
  s.l = 1,
  s.b = 1,
  rho = c(0.5, 0.5, 0.5),
  beta = rep(0, length(beta.truth))
)
tuning_parameters <- list(
  s.l.tune = 0.1,
  s.b.tune = 0.1,
  rho.tune = c(0.1, 0.1, 0.1)
)
prior_parameters <- list(
  s.l.prior = c(1, 1),
  s.b.prior = c(1, 1),
  rho.prior = c(1, 1, 1),
  mu.beta = matrix(0, length(beta.truth), 1),
  s2.beta = diag(length(beta.truth))
)

# Source the MCMC_Algorithm.R script from GitHub
source(paste0('https://raw.githubusercontent.com/perrywilliamsunr/',
              'Wood-Duck-Mercury-Contamination/main/Code/MCMC_Algorithm/',
              'MCMC_Algorithm.R'))

# Run the MCMC sampling
run_mcmc(y, X, n.iter, initial_values, tuning_parameters, prior_parameters, file_name, save_dir)
# After running this script, you should find output files in the specified directory

# Load the MCMC output file
load(file.path(save_dir, paste0(file_name, ".iter", n.iter, ".RData")))

# Plotting MCMC results
par(mfrow=c(2, 2))
for (i in 1:ncol(MCMC.Output$rho)) {
  plot(MCMC.Output$rho[, i], type = 'l', ylim = c(0, 1), main = paste("Trace plot of rho[", i, "]", sep = ""))
}
q <- apply(MCMC.Output$rho, 2, quantile, c(0.025, 0.975))

par(mfrow = c(3, 3))
for (i in 1:ncol(MCMC.Output$beta)) {
  plot(MCMC.Output$beta[, i], type = 'l', ylim = c(-10, 10), main = paste("Trace plot of beta[", i, "]", sep = ""))
  # Add lines for true values
  abline(h = beta.truth[i], col = 2)
}
q <- apply(MCMC.Output$beta, 2, quantile, c(0.025, 0.975))
q[1, ] < beta.truth & q[2, ] > beta.truth

par(mfrow = c(2, 1))
plot(MCMC.Output$s.l, type = 'l', main = "Trace plot of s.l")
# Add line for true value if available
abline(h = s.l.truth, col = 2)

par(mfrow = c(2, 1))
plot(MCMC.Output$s.b, type = 'l', main = "Trace plot of s.b")
# Add line for true value if available
abline(h = s.b.truth, col = 2)

