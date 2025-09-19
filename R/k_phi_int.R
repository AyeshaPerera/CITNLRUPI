k_phi_int <- function(rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delt, no.nodes.GL){
  # This module computes
  #        delta        d
  #    integral   integral  k(h, gamma, rho)
  #     -delta       -d
  #          * phi(h - gamma) dh{p-1} dh{p}
  # using Gauss quadrature
  #
  # Inputs
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
  # delt: tuning constant for w1(x) were 0 < delt <= d.
  # no.nodes.GL: the number of Gauss Legendre nodes
  #
  # Output
  #        delta        d
  #    integral   integral  k(h, gamma, rho)
  #     -delta       -d
  #         * phi(h - gamma) dh{p-1} dh{p}
  #
  # Written by A.Perera

  GL.modified.2 <- GL_modified_nodes_weights(-d, d, no.nodes.GL)
  GL.modified.1 <- GL_modified_nodes_weights(-delt, delt, no.nodes.GL)
  cp <- 0

  for(j in 1:no.nodes.GL){
    h.p <- GL.modified.1$GL.modified.nodes[j]
    w.h.p <- GL.modified.1$GL.modified.weights[j]

    for (i in 1:no.nodes.GL) {
      h.p.1 <- GL.modified.2$GL.modified.nodes[i]
      w.h.p.1 <- GL.modified.2$GL.modified.weights[i]
      k.x <- k_fn(h.p.1,  h.p, rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delt)
      cp.val <- (k.x) * dnorm(h.p.1 - gamma.vec[1]) * dnorm(h.p - gamma.vec[2])
      cp <- cp + w.h.p.1 * w.h.p * cp.val
    }
  }
  cp
}
