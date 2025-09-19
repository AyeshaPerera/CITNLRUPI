CP_MC_vr <- function(gamma.vec, rho, N, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta){
  
  set.seed(890)
  # This module calculates the coverage probability using Monte  
  # Carlo simulation variance reduction by conditioning
  #
  # Input
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
  
  # Simulate random vector [H] using multivariate normal distribution which (gamma, I)
  obs.h <- rmvnorm(N, mean = gamma.vec, sigma = matrix(c(1, 0, 0, 1), nrow = 2))
  
  k.obs.h <- rep(0, N)
  for (i in 1:N) {
    k.obs.h[i] <- k_fn(obs.h[i, 1], obs.h[i, 2], rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta)
  }
  
  cp.estimate <- mean(k.obs.h)
  s.sq <- sum((k.obs.h - cp.estimate)^2) / (N - 1)
  se <- sqrt(s.sq / N)
  
  # List to store the estimated coverage probability and it's standard error
  list(Estimate=cp.estimate, SE=se)
  
}
