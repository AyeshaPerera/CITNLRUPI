SEL_i <- function(r.knots.vec, s.knots.vec, t.vec, alpha, d, gamma.vec, delta){
  # This module computes the scaled expected length
  # using integrate()
  #
  # Inputs
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # no.nodes.GL: the number of Gauss Legendre nodes
  #
  # Output
  # Scaled expected length using integrate()
  #
  # Written by A.Perera August 2023
  
  1 + ISEL1_i(s.knots.vec, t.vec, alpha, d, gamma.vec) + ISEL2_i(r.knots.vec, t.vec, alpha, d, gamma.vec, delta) - (pnorm(d - gamma.vec[2]) - pnorm(-d - gamma.vec[2]))
}
