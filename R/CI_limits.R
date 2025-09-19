CI_limits <- function(h, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta){
  
  # This module computes the endpoints of the confidence 
  # interval where
  # l(h) = b(hp) - s(hp) - w1(|hp|) * (b1(h{p-1}) + r(h{p-1}))
  # u(h) = b(hp) + s(hp) - w1(|hp|) * (b1(h{p-1}) - r(h{p-1}))
  #
  # Inputs
  # h: the vector (hp, h{p-1})
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # delta: tuning constant for w1(x) were 0 < delt <= d. 
  #
  # Output
  # upper and lower limits
  #
  # Written by A.Perera May 2023 
  
  b.hp <- b_fn(h[2], b.knots.vec, t.vec, d)
  w1.hp <- w1_fn(abs(h[2]), d, delta)
  b1.hp1 <- b1_fn(h[1], b1.knots.vec, t.vec, d)
  s.hp <- s_fn(h[2], s.knots.vec, t.vec, alpha, d)
  r.hp1 <- r1_fn(h[1], r.knots.vec, t.vec, alpha, d)
  
  uh <- b.hp + s.hp + w1.hp * (b1.hp1 + r.hp1)
  lh <- b.hp - s.hp + w1.hp * (b1.hp1 - r.hp1)
  
  limits <- list("l" = lh, "u" = uh)  
  limits
}
