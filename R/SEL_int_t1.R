SEL_int_t1 <- function(s.knots.vec, t.vec, alpha, d, gamma.vec){
  # This module computes
  #     1           d
  # ------------ integrate s(h{p}) * psi(h{p} - gamma{p}) dh{p}
  # z(1-alpha/2)   -d
  # using the integrate function
  #
  # Inputs
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  #
  # Output
  # first term integration part of the equation to scaled expected length
  #
  # Written by A.Perera May 2023
  
  z <- qnorm(1-alpha/2)
  
  fun <- function(h.p){
    # This module represents s(h{p}) * psi(h{p} - gamma{p})
    
    z <- qnorm(1-alpha/2)
    s.h.p <- s_fn(h.p, s.knots.vec, t.vec, alpha, d)
    (s.h.p)* dnorm(h.p - gamma.vec[2])
  }
  
  isel <- integrate(function(h.p){sapply(h.p, fun)}, lower = -d, upper = d, rel.tol = 1e-10)$value
  
  isel / z
  
}