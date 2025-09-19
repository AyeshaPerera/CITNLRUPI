SEL_int_t2 <- function(r.knots.vec, t.vec, alpha, d, gamma.vec, delta){
  # This module computes
  #     1           d
  # ------------ (integrate (w1(|h{p}|))psi(h{p} - gamma{p}) dh{p}
  # z(1-alpha/2)   -d
  #     d
  # integrate r1(h{p-1}) psi(h{p-1} - gamma{p-1}) dh{p-1})
  #   -d
  # using the integrate function
  #
  # Inputs
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # delta: 
  #
  # Output
  # second term integration part of the equation to scaled expected length
  #
  # Written by A.Perera May 2023
  
  z <- qnorm(1-alpha/2)
  
  fun1 <- function(h.p){
    # This module represents (w1(|h{p}|))psi(h{p} - gamma{p})
    w1_fn(abs(h.p), d, delta) * dnorm(h.p - gamma.vec[2])
  }
  
  isel1 <- integrate(function(h.p){sapply(h.p, fun1)}, -delta, delta)$value
  
  fun2 <- function(h.p.1){
    # This module represents r1(h{p-1}) psi(h{p-1} - gamma{p-1})
    r1_fn(h.p.1, r.knots.vec, t.vec, alpha, d) * dnorm(h.p.1 - gamma.vec[1])
  }
  
  isel2 <- integrate(function(h.p.1){sapply(h.p.1, fun2)}, lower = -d, upper = d, rel.tol = 1e-10)$value
  isel1 * isel2/ z
  
}
