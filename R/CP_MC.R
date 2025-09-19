library(mvtnorm)
CP_MC <- function(gamma.vec, rho, N, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta) {
  
  # This module calculates the coverage probability using 
  # Monte Carlo simulation
  #
  # Inputs
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # rho: vector of (rho{p-1}, rho{p})
  # N: Number of independent simulation runs
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # delta: tuning constant for w1(x) were 0 < delta <= d
  #
  # Output:
  # coverage probability estimate
  #
  # Written by A.Perera June 2023
  
  # Define parameters for multivariate normal distribution
  mean.vec <- c(0, gamma.vec)
  cov.mat <- matrix(c(1, rho[1], rho[2], rho[1], 1, 0, rho[2], 0, 1), nrow=3)
  
  # Simulate random vector [G, H] using multivariate normal distribution
  sim.data <- rmvnorm(N, mean = mean.vec, sigma = cov.mat)
  
  # Extract G and H from simulated data
  G <- sim.data[, 1]
  H <- sim.data[, -1]
  
  # Calculate indicator function for each simulation run
  indicators <- rep(0, N)
  for(i in 1:N){
    limits.list <- CI_limits(H[i,], b.knots.vec, b1.knots.vec, 
                             s.knots.vec, r.knots.vec, t.vec, alpha, d, delta)
    if(limits.list$l <= G[i] & G[i] <= limits.list$u){
      indicators[i] <- 1
    }else{
      indicators[i] <- 0
    }
    
  }
  
  # Estimate coverage probability using Monte Carlo simulation
  cp.estimate <-  sum(indicators) / N
  
  # Standard error
  se <- sqrt(cp.estimate * (1 - cp.estimate) / N)
  
  # List to store the estimated coverage probability and it's standard error
  list(Estimate=cp.estimate, SE=se)
}
