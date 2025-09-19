SEL_MC_vr <- function(gamma.vec, N, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta){
  
  set.seed(100)
  
  # This module calculates the scaled expected length using Monte Carlo simulation with common random numbers
  #
  # Input
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # N: Number of independent simulation runs
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # delta: tuning constant for w1(x) were 0 < delta <= d 
  #
  # Output:
  # Scaled expected length estimate
  #
  # Written by A.Perera June 2023
  
  # List to store the estimated scaled expected length and it's standard error
  sel.list <- list()
  
  # generate random values
  values <- rmvnorm(N*length(gamma.vec), mean = gamma.vec, sigma = matrix(c(1, 0, 0, 1), nrow = 2))
  
  func <- function(h, s.knots.vec, r.knots.vec, t.vec, alpha, d){
    
    z <- qnorm(1 - alpha/2)
    h.p <- h[2]
    h.p.1 <- h[1]
    s.h.p <- s_x(h.p, s.knots.vec, t.vec, alpha, d) 
    w.1.h.p <- w_1_x(abs(h.p), d, delta) 
    r.1.h.p <- r_1_x(h.p.1, r.knots.vec, t.vec, alpha, d)
    (s.h.p + w.1.h.p * r.1.h.p) / z
  }
  
  estimates <- apply(values, 1, func, s.knots.vec=s.knots.vec, r.knots.vec=r.knots.vec, t.vec=t.vec, alpha=alpha, d=d)
  mc.estimate <- mean(estimates)
  s <- sqrt(sum((estimates - mc.estimate)^2) / (N - 1))
  se <- s / sqrt(N)
  sel.list$Estimate <- mc.estimate
  sel.list$SE <- se
  sel.list
}
