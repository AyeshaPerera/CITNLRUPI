ICP_T1_GL <- function(rho, d, alpha, gamma.vec, t.vec, 
         b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, 
         delta, no.nodes.GL){  
  # This module computes 
  #    d        Inf
  #  integral integral  (k(h, gamma, rho) - k+(h, gamma, rho)) 
  #     -d       -Inf
  # 	* phi(h - gamma) dh{p-1} dh{p}
  # using Gauss Legendre quadrature 
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
  # 	* phi(h - gamma) dh{p-1} dh{p}
  #
  # Written by A.Perera
  
  GL.modified <- GL_modified_nodes_weights(-d, d,
           no.nodes.GL)
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
  
  GL.modified.2.1 <- GL_modified_nodes_weights(0, 1, 
           no.nodes.GL)
  term.2.1 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    
    for(i in 1:no.nodes.GL){
      t <- GL.modified.2.1$GL.modified.nodes[i]
      h.p.1.i <- d + (1 - t)/t
      w.h.p.1.i <- GL.modified.2.1$GL.modified.weights[i]
      term1.val <- (1/(t^2))*I_fn(h.p.1.i, h.p.j, rho, d, 
               alpha, gamma.vec, t.vec, b.knots.vec, 
               b1.knots.vec, s.knots.vec, r.knots.vec, delta)
      term.2.1 <- term.2.1 + w.h.p.j * w.h.p.1.i * term1.val
    }
  }
  
  GL.modified.2.2 <- GL_modified_nodes_weights(0, 1, 
           no.nodes.GL)
  term.2.2 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    
    for(i in 1:no.nodes.GL){
      t <- GL.modified.2.2$GL.modified.nodes[i]
      h.p.1.i <- -d - (1 - t)/t
      w.h.p.1.i <- GL.modified.2.2$GL.modified.weights[i]
      term2.val <- (1/(t^2))*I_fn(h.p.1.i, h.p.j, rho, d, alpha, 
               gamma.vec, t.vec, b.knots.vec, b1.knots.vec, 
               s.knots.vec, r.knots.vec, delta)
      term.2.2 <- term.2.2 + w.h.p.j * w.h.p.1.i * term2.val
    }
  }
  term.1 + term.2.1 + term.2.2
}
