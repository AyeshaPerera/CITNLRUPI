ICP_T1_int <- function(rho, d, alpha, gamma.vec, t.vec, 
         b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, 
         delta, no.nodes.GL){
  
  # This module computes 
  #    d        Inf
  #  integral integral  (k(h, gamma, rho) - k+(h, gamma, rho))
  #     -d       -Inf
  #         * phi(h - gamma) dh{p-1} dh{p}
  # using integrate function
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
  # delta: tuning constant for w1(x) were 0 < delta <= d
  #
  # Output
  #    d        Inf
  #  integral integral  (k(h, gamma, rho) - k+(h, gamma, rho))
  #    -d       -Inf
  #         * phi(h - gamma) dh{p-1} dh{p}
  #
  # Written by A.Perera
  
  GL.modified <- GL_modified_nodes_weights(-d, d, no.nodes.GL)
  term.1 <- 0
  
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    for (i in 1:no.nodes.GL) {
      h.p.1.i <- GL.modified$GL.modified.nodes[i]
      w.h.p.1.i <- GL.modified$GL.modified.weights[i]
      term.val <- I_fn(h.p.1.i, h.p.j, rho, d, alpha, gamma.vec, 
               t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, 
               r.knots.vec, delta)
      term.1 <- term.1 + w.h.p.1.i * w.h.p.j * term.val
      
    }
  }
  
  integral1 <- function(h.p) {
    integrate(function(h.p.1){sapply(h.p.1, I_fn,h.p=h.p, 
             rho=rho, d=d, alpha=alpha, gamma.vec=gamma.vec, 
             t.vec=t.vec, delta=delta, b.knots.vec=b.knots.vec, 
             b1.knots.vec=b1.knots.vec, s.knots.vec=s.knots.vec,
             r.knots.vec=r.knots.vec)},lower = -Inf, upper = -d, 
             rel.tol = 1e-15)$value 
  }
  
  term.2.3.1 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.3.1 <- term.2.3.1 + w.h.p.j * integral1(h.p.j)
    
  }
  
  integral2 <- function(h.p) {
    integrate(function(h.p.1){sapply(h.p.1, I_fn,h.p=h.p, 
             rho=rho, d=d, delta=delta, alpha=alpha, 
             gamma.vec=gamma.vec, t.vec=t.vec, 
             b.knots.vec=b.knots.vec, b1.knots.vec=b1.knots.vec, 
             s.knots.vec=s.knots.vec, r.knots.vec=r.knots.vec)}, 
             lower = d, upper = Inf, rel.tol = 1e-15)$value 
  }
  
  term.2.3.2 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.3.2 <- term.2.3.2 + w.h.p.j * integral2(h.p.j)
  }
  
  term.1 + term.2.3.1 + term.2.3.2
  
}
