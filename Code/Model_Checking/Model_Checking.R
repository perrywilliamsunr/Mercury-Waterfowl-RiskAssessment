# Clear the workspace
rm(list=ls())

#######################################################
### Script to examine Bayesian p-values of each model
#######################################################

# Load necessary libraries
library(coda)        # For trace plots and density plots
library(truncnorm)   # For truncnorm function
library(MASS)        # For mvrnorm function
library(Matrix)      # For Kronecker function
library(mvtnorm)     # For multivariate normal distribution
library(R.utils)     # For saving files

# Define your data and parameters
n.iter <- 10000  # Number of MCMC iterations
file_name <- "Bayes_PValue_Output"  # Base name for saving files
save_dir <- "/Users/pwill/Dropbox/GitHub/Wood-Duck-Mercury-Contamination/R_Output"  # Directory to save output

# Load data from GitHub
data <- read.csv(paste0("https://github.com/perrywilliamsunr/",
                        "Wood-Duck-Mercury-Contamination/raw/main/Data/data.csv"))

# Prepare data for MCMC
n <- length(unique(data$bird_id))
ind <- seq(1, n * 3, 3)
y <- c(data$y1[ind[1]:(ind[1] + 2)], data$y2[ind[1]:(ind[1] + 2)])
for(i in 2:length(ind)){
  y.tmp <- c(data$y1[ind[i]:(ind[i] + 2)], data$y2[ind[i]:(ind[i] + 2)])
  y <- c(y, y.tmp)
}
y.untransformed <- y
y <- log(y)  # Log-transform the data

# For reproducibility
set.seed(123)

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
m.tmp5[c(4, 5, 6, 10, 11, 12, 16, 17, 18)] <- 1
m.tmp5[c(19, 20, 21, 25, 26, 27, 31, 32, 33)] <- 1

# Independent variables
n.mall <- 5
liver <- rep(rep(1:0, each = 3), n)
indices <- seq(1, nrow(data), by = 3)
species <- rep(data$species, each = 2)
flank.tmp <- (data$flank1 + data$flank2) / 2
flank.scaled <- (flank.tmp - 8821.888) / 9209.355  # Scaling flank
axillary.tmp <- (data$axillary1 + data$axillary2) / 2
axillary.scaled <- (axillary.tmp - mean(unique(axillary.tmp))) / sd(unique(axillary.tmp))  # Scaling axillary
breast.tmp <- (data$breast1 + data$breast2) / 2
breast.scaled <- (breast.tmp - mean(unique(breast.tmp))) / sd(unique(breast.tmp))  # Scaling breast
flank <- rep(flank.scaled, each = 2)
axillary <- rep(axillary.scaled, each = 2)
breast <- rep(breast.scaled, each = 2)

# Model selection variables
model.number <- 35
Species <- TRUE
Liver <- TRUE
Breast <- FALSE
Flank <- TRUE
Species.Liver <- TRUE
Species.Flank <- TRUE
Species.Axillary <- FALSE
Species.Breast <- FALSE

# Construct design matrix X
X <- matrix(1, length(breast), 1)
if(Liver) X <- cbind(X, liver)
if(Flank) X <- cbind(X, flank)
if(Breast) X <- cbind(X, breast)
if(Species) X <- cbind(X, species)
if(Species.Liver) X <- cbind(X, species * liver)
if(Species.Flank) X <- cbind(X, species * flank)
if(Species.Axillary) X <- cbind(X, (1 - species) * axillary)
if(Species.Breast) X <- cbind(X, species * breast)

# Set initial values, tuning parameters, and prior parameters
initial_values <- list(
  s.l = 1,
  s.b = 1,
  rho = c(0.5, 0.5, 0.5),
  beta = rep(0, ncol(X))
)
tuning_parameters <- list(
  s.l.tune = 0.1,
  s.b.tune = 0.1,
  rho.tune = c(0.1, 0.1, 0.1)
)
prior_parameters <- list(
  s.l.prior = c(825, 100^2),
  s.b.prior = c(.01, .01),
  rho.prior = c(100, 1),
  mu.beta = matrix(0, ncol(X), 1),
  s2.beta = 1 * diag(ncol(X))
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
burn <- floor(n.iter / 2)
thin <- 1
ind <- seq(burn, n.iter, thin)
y.ppd.save <- MCMC.Output$y.ppd[ind,]
beta.save <- MCMC.Output$beta[ind,]
rho.save <- MCMC.Output$rho[ind,]
s.l.save <- MCMC.Output$s.l[ind,]
s.b.save <- MCMC.Output$s.b[ind,]
Dbar.save <- MCMC.Output$Dbar[ind,]

# Initialize matrices to store Bayesian p-value counts
Bayes.count1 <- matrix(0, length(ind), length(y))
Bayes.count2 <- matrix(0, length(ind), 1)

# Calculate Bayesian p-values
for(k in 1:length(ind)){
  mu <- matrix(X %*% beta.save[k,], n * 6, 1)
  block <- s.l.save[k]^2 * m.tmp1 +
    s.b.save[k]^2 * m.tmp2 +
    rho.save[k, 1] * m.tmp3 * s.l.save[k]^2 +
    rho.save[k, 2] * m.tmp4 * s.b.save[k]^2 +
    rho.save[k, 3] * s.l.save[k] * s.b.save[k] * m.tmp5
  Sigma <- kronecker(diag(n), block)
  Bayes.count1[k, ] <- ifelse((log(y.ppd.save[k,]) - mu)^2 / mu < 
                                ((y - mu)^2 / mu), 1, 0)
  if(-2 * dmvnorm(log(y.ppd.save[k,]), mu, Sigma, log = TRUE) < 
     -2 * dmvnorm(y, mu, Sigma, log = TRUE)){
    Bayes.count2[k, 1] <- 1
  }
}

# Calculate and print Bayesian p-values
Bayes.p <- mean(Bayes.count2)
print(Bayes.p)

# Plot trace plots for parameters
par(mfrow = c(2, 2))
for(i in 1:ncol(MCMC.Output$rho)){
  plot(MCMC.Output$rho[, i], type = 'l', ylim = c(0, 1), main = paste("Trace plot of rho[", i, "]", sep = ""))
}

par(mfrow = c(3, 3))
for(i in 1:ncol(MCMC.Output$beta)){
  plot(MCMC.Output$beta[, i], type = 'l', ylim = c(-10, 10), main = paste("Trace plot of beta[", i, "]", sep = ""))
}

par(mfrow = c(2, 1))
plot(MCMC.Output$s.l, type = 'l', main = "Trace plot of s.l")

par(mfrow = c(2, 1))
plot(MCMC.Output$s.b, type = 'l', main = "Trace plot of s.b")
