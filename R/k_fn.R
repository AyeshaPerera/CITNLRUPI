k_fn <- function(h.p.1, h.p, rho, d, alpha, gamma.vec, t.vec, 
         b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta){
  # This module computes k(h, gamma, rho) where
  # k(x, w, gamma, rho) = Psi(l(h), u(h); 
  #         rho(h - gamma), 1 - rho^2)
  # 
  # Inputs
  # h.p.1: a value (h{p-1} of h vector)
  # h.p: a value (h{p} of h vector)
  # rho: vector of (rho{p-1}, rho{p})
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # delta: tuning constant for w1(x) were 0 < delta <= d. 
  #
  # Output
  # Value of k(h, gamma, rho)
  #
  # Written by A.Perera May 2023
  
  h <- c(h.p.1, h.p)
  limits <- CI_limits(h, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta)
  rho.col <- matrix(rho)
  h.col <- matrix(h)
  gamma.vec.col <- matrix(gamma.vec)
  mu <- t(rho.col) %*% (h.col - gamma.vec.col)
  variance <- 1 - t(rho.col) %*% rho.col
  
  Psi(limits$l, limits$u, mu, variance)
}
