run_mcmc <- function(y, X, n.iter, initial_values, tuning_parameters, prior_parameters, file_name, save_dir) {
  
  n=nrow(X)/6
  
  # Extract initial values
  s.l <- initial_values$s.l
  s.b <- initial_values$s.b
  rho <- initial_values$rho
  beta <- initial_values$beta
  
  # Extract tuning parameters
  s.l.tune <- tuning_parameters$s.l.tune
  s.b.tune <- tuning_parameters$s.b.tune
  rho.tune <- tuning_parameters$rho.tune
  
  # Extract prior parameters
  s.l.prior <- prior_parameters$s.l.prior
  s.b.prior <- prior_parameters$s.b.prior
  rho.prior <- prior_parameters$rho.prior
  mu.beta <- prior_parameters$mu.beta
  s2.beta <- prior_parameters$s2.beta
  
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
  
  # Starting values
  s.l <- 1
  s.b <- 1
  rho <- rep(0.5, 3)
  beta <- matrix(0, ncol(X), 1)
  
  beta <- matrix(beta, length(beta), 1)
  mu <- matrix(X %*% beta, n * 6, 1)
  block <- s.l^2 * m.tmp1 +
    s.b^2 * m.tmp2 +
    rho[1] * m.tmp3 * s.l^2 +
    rho[2] * m.tmp4 * s.b^2 +
    rho[3] * s.l * s.b * m.tmp5
  Sigma <- kronecker(diag(n), block)
  
  # Book keeping
  accept.s.l <- 0
  accept.s.b <- 0
  accept.rho <- rep(0, 3)
  s.l.save <- matrix(NA, n.iter, 1)
  s.l.tune.save <- s.l.tune
  s.b.save <- matrix(NA, n.iter, 1)
  s.b.tune.save <- s.b.tune
  rho.save <- matrix(NA, n.iter, 3)
  beta.save <- matrix(NA, n.iter, length(beta))
  rho.tune.save <- matrix(NA, n.iter, 3)
  y.ppd.save <- matrix(NA, n.iter, length(y))
  Dbar.save <- matrix(NA, n.iter, 1)
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  
  # MCMC loop
  for (k in 1:n.iter) {
    # Update progress bar
    setTxtProgressBar(pb, k)
    
    # Sample s.l
    s.l.star <- rtruncnorm(1, a = 0, b = 10000, s.l, s.l.tune)
    block <- s.l.star^2 * m.tmp1 +
      s.b^2 * m.tmp2 +
      rho[1] * m.tmp3 * s.l.star^2 +
      rho[2] * m.tmp4 * s.b^2 +
      rho[3] * s.l.star * s.b * m.tmp5
    Sigma.star <- kronecker(diag(n), block)
    mh1 <- dmvnorm(c(y), c(mu), Sigma.star, log = TRUE) +
      dinvgamma(s.l.star, shape = s.l.prior[1], rate = s.l.prior[2], log = TRUE) +
      log(dtruncnorm(s.l, a = 0, b = 10000, mean = s.l.star, sd = s.l.tune))
    mh2 <- dmvnorm(c(y), c(mu), Sigma, log = TRUE) +
      dinvgamma(s.l, shape = s.l.prior[1], rate = s.l.prior[2], log = TRUE) +
      log(dtruncnorm(s.l.star, a = 0, b = 10000, mean = s.l, sd = s.l.tune))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      s.l <- s.l.star
      Sigma <- Sigma.star
      accept.s.l <- accept.s.l + 1
    }
    
    # Sample s.b
    s.b.star <- rtruncnorm(1, a = 0, b = 10000, s.b, s.b.tune)
    block <- s.l^2 * m.tmp1 +
      s.b.star^2 * m.tmp2 +
      rho[1] * m.tmp3 * s.l^2 +
      rho[2] * m.tmp4 * s.b.star^2 +
      rho[3] * s.l * s.b.star * m.tmp5
    Sigma.star <- kronecker(diag(n), block)
    mh1 <- dmvnorm(c(y), c(mu), Sigma.star, log = TRUE) +
      dinvgamma(s.b.star, shape = s.b.prior[1], rate = s.b.prior[2], log = TRUE) +
      log(dtruncnorm(s.b, a = 0, b = 10000, mean = s.b.star, sd = s.b.tune))
    mh2 <- dmvnorm(c(y), c(mu), Sigma, log = TRUE) +
      dinvgamma(s.b, shape = s.b.prior[1], rate = s.b.prior[2], log = TRUE) +
      log(dtruncnorm(s.b.star, a = 0, b = 10000, mean = s.b, sd = s.b.tune))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      s.b <- s.b.star
      Sigma <- Sigma.star
      accept.s.b <- accept.s.b + 1
    }
    
    # Sample rho1
    rho1.star <- rtruncnorm(1, a = 0, b = 1, rho[1], rho.tune[1])
    block <- s.l^2 * m.tmp1 +
      s.b^2 * m.tmp2 +
      rho1.star * m.tmp3 * s.l^2 +
      rho[2] * m.tmp4 * s.b^2 +
      rho[3] * s.l * s.b * m.tmp5
    Sigma.star <- kronecker(diag(n), block)
    mh1 <- dmvnorm(c(y), c(mu), Sigma.star, log = TRUE) +
      sum(log(dtruncnorm(rho[1], a = 0, b = 1, mean = rho1.star, sd = rho.tune[1])))
    mh2 <- dmvnorm(c(y), c(mu), Sigma, log = TRUE) +
      sum(log(dtruncnorm(rho1.star, a = 0, b = 1, mean = rho[1], sd = rho.tune[1])))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      rho[1] <- rho1.star
      Sigma <- Sigma.star
      accept.rho[1] <- accept.rho[1] + 1
    }
    
    # Sample rho2
    rho2.star <- rtruncnorm(1, a = 0, b = 1, rho[2], rho.tune[2])
    block <- s.l^2 * m.tmp1 +
      s.b^2 * m.tmp2 +
      rho[1] * m.tmp3 * s.l^2 +
      rho2.star * m.tmp4 * s.b^2 +
      rho[3] * s.l * s.b * m.tmp5
    Sigma.star <- kronecker(diag(n), block)
    mh1 <- dmvnorm(c(y), c(mu), Sigma.star, log = TRUE) +
      sum(log(dtruncnorm(rho[2], a = 0, b = 1, mean = rho2.star, sd = rho.tune[2])))
    mh2 <- dmvnorm(c(y), c(mu), Sigma, log = TRUE) +
      sum(log(dtruncnorm(rho2.star, a = 0, b = 1, mean = rho[2], sd = rho.tune[2])))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      rho[2] <- rho2.star
      Sigma <- Sigma.star
      accept.rho[2] <- accept.rho[2] + 1
    }
    
    # Sample rho3
    rho3.star <- rtruncnorm(1, a = 0, b = 1, rho[3], rho.tune[3])
    block <- s.l^2 * m.tmp1 +
      s.b^2 * m.tmp2 +
      rho[1] * m.tmp3 * s.l^2 +
      rho[2] * m.tmp4 * s.b^2 +
      rho3.star * s.l * s.b * m.tmp5
    Sigma.star <- kronecker(diag(n), block)
    mh1 <- dmvnorm(c(y), c(mu), Sigma.star, log = TRUE) +
      sum(log(dtruncnorm(rho[3], a = 0, b = 1, mean = rho3.star, sd = rho.tune[3])))
    mh2 <- dmvnorm(c(y), c(mu), Sigma, log = TRUE) +
      sum(log(dtruncnorm(rho3.star, a = 0, b = 1, mean = rho[3], sd = rho.tune[3])))
    mh <- exp(mh1 - mh2)
    if (mh > runif(1)) {
      rho[3] <- rho3.star
      Sigma <- Sigma.star
      accept.rho[3] <- accept.rho[3] + 1
    }
    
    # Sample beta (conjugate)
    A <- t(X) %*% solve(Sigma) %*% X + solve(s2.beta)
    b.t <- t(y) %*% solve(Sigma) %*% X + t(mu.beta) %*% solve(s2.beta)
    beta <- rmvnorm(1, solve(A) %*% t(b.t), solve(A))
    
    # Model checking
    mu <- matrix(X %*% t(beta), n * 6, 1)
    y.ppd.save[k,] <- matrix(exp(rmvnorm(1, mean = mu, Sigma)), , 1)
    
    # DIC calculations
    Dbar.save[k, 1] <- -2 * dmvnorm(c(y), c(X %*% t(beta)), Sigma, log = TRUE)
    
    # Auto tune
    if (accept.s.l / k < 0.3) s.l.tune <- s.l.tune * 0.9
    if (accept.s.l / k > 0.5) s.l.tune <- s.l.tune * 1.1
    s.l.tune.save[k] <- s.l.tune
    
    if (accept.s.b / k < 0.3) s.b.tune <- s.b.tune * 0.9
    if (accept.s.b / k > 0.5) s.b.tune <- s.b.tune * 1.1
    s.b.tune.save[k] <- s.b.tune
    
    rho.tune <- ifelse(accept.rho / k < 0.3, rho.tune * 0.9, rho.tune)
    rho.tune <- ifelse(accept.rho / k > 0.5, rho.tune * 1.1, rho.tune)
    rho.tune.save[k,] <- rho.tune
    
    # Store samples
    beta.save[k,] <- beta
    s.l.save[k] <- s.l
    s.b.save[k] <- s.b
    rho.save[k,] <- rho
    
    # Save results periodically
    if (k %% (n.iter * 0.25) == 0) {
      MCMC.Output <- list(
        beta = beta.save,
        rho = rho.save,
        s.l = s.l.save,
        s.b = s.b.save,
        Dbar = Dbar.save,
        y.ppd = y.ppd.save
      )
      
      # Create directory if it doesn't exist
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
      }
      
      # Save file
      save(MCMC.Output, file = paste0(save_dir, "/", file_name, ".iter", k, ".RData"))
    }
  }
  
  # Close progress bar
  close(pb)
}
