CP_T1_mult_int <- function(gamma.vec, alpha, d, rho, t.vec,
         s.knots.vec, b.knots.vec, b1.knots.vec, r.knots.vec, 
         delta, no.nodes.GL){
  
  # This module computes the coverage probability
  #
  # Inputs
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # rho: vector of (rho{p-2}, rho{p})
  # t.vec: the vector (t0, t1, ..., tk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  #
  # Output
  # Coverage probability computed with integrate functions
  #
  # Written by A.Perera
  
  cp <- (1 - alpha) + ICP_T1_mult_int(rho, d, alpha, gamma.vec, 
           t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, 
           r.knots.vec, delta, no.nodes.GL)
  
  cp
  
}