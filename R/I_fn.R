I_fn <- function(h.p.1, h.p, rho, d, alpha, gamma.vec, t.vec,
         b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta){
  # This module represents
  # (k(h, gamma, rho) - k+(h, gamma, rho)) * phi(h - gamma)
  # using integrate function
  #
  # Inputs
  # h.p.1: a value
  # h.p: a value
  # rho: vector of (rho{p-1}, rho{p})
  # d: a positive number that specifies
  #    the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # delta: tuning constant for w1(x) were 0 < delta <= d.
  #
  # Output
  # represents (k(h, gamma, rho) - k+(h, gamma, rho))
  #         * phi(h - gamma)
  #
  # Written by A.Perera

  k.x <- k_fn(h.p.1,  h.p, rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta)
  k.dag.x <- k_dag_fn(h.p.1, h.p, rho, d, alpha, gamma.vec)
  (k.x - k.dag.x) * dnorm(h.p.1 - gamma.vec[1]) * dnorm(h.p - gamma.vec[2])
}
