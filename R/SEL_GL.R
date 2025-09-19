SEL_GL <- function(r.knots.vec, s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta){
  # This module computes the scaled expected length
  # using Gauss Legendre quadrature
  #
  # Inputs
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # no.nodes.GL: the number of Gauss Legendre nodes
  #
  # Output
  # Scaled expected length
  #
  # Written by A.Perera August 2023
  
  1 + SEL_GL_t1(s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL) + SEL_GL_t2(r.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta) - (pnorm(d - gamma.vec[2]) - pnorm(-d - gamma.vec[2]))
}