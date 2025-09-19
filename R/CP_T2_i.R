CP_T2_int <- function(rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta){
  # Computes coverage probability using integrate function 
  #
  # Input
  # rho: vector of (rho{p-2}, rho{p})
  # d: a positive number that specifies the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # delta: tuning constant for w1(x) were 0 < delta <= d 
  # no.nodes.GL
  #
  # Output
  # Coverage probability
  #
  # Written by A. Perera June 2023
  
  z <- qnorm(1 - alpha / 2)
  mu <- c(0, gamma.vec[2])
  Sigma <- matrix(c(1, rho[2], rho[2], 1), nrow=2)
  
  term1 <- pmvnorm(lower=c(-z, -d), upper=c(z, d), 
           algorithm = "Miwa", mean=mu, sigma=Sigma)[1]
  
  term2 <- k_phi_int_i(rho, d, alpha, gamma.vec, t.vec, 
           b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, 
           delta)
  
  term3 <- CP_T2_int_t3(rho, gamma.vec, alpha, d, t.vec, 
           b.knots.vec, s.knots.vec, delta)
  
  term4 <- CP_T2_int_t4(rho, d, alpha, gamma.vec, t.vec, 
           b.knots.vec, s.knots.vec, delta)
  
  1 - alpha - term1 + term2 + term3 + term4
  
}
