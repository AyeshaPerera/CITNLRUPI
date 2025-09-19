CP_T2_GL_t4 <- function(rho, gamma.p, alpha, d, t.vec, 
         b.knots.vec, s.knots.vec, delt, no.nodes.GL){
  # This module computes
  #                      -delt    1
  #  t_4(gamma_p; b, s) = int   int  (f_{G, H_p}(b(h_p) 
  #                      -d     -1
  #      + s(h_p)y, h_p) * s(h_p) dy dh_p)
  #      d     1
  #  +  int   int  (f_{G, H_p}(b(h_p) + s(h_p)y, h_p)
  #   delt   -1
  #       * s(h_p) dy dh_p)
  #
  # Input
  # rho: vector of (rho{p-2}, rho{p})
  # d: a positive number that specifies the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.p: gamma{p} value
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # delt: tuning constant for w1(x) were 0 < delt <= d 
  #
  # Output
  #                      -delt    1
  #  t_4(gamma_p; b, s) = int   int  (f_{G, H_p}(b(h_p) 
  #                      -d     -1
  #      + s(h_p)y, h_p) * s(h_p) dy dh_p)
  #      d     1
  #  +  int   int  (f_{G, H_p}(b(h_p) + s(h_p)y, h_p)
  #   delt   -1
  #       * s(h_p) dy dh_p)
  #
  # Written by A.Perera June 2023
  
  fn <- function(h.p, y, rho, gamma.p, alpha, d, t.vec, 
           b.knots.vec, s.knots.vec){
    b.h.p <- b_fn(t=h.p, b.knots.vec=b.knots.vec, t.vec=t.vec, 
             d=d)
    s.h.p <- s_fn(t=h.p, s.knots.vec=s.knots.vec, t.vec=t.vec, 
             alpha = alpha, d=d)
    g <- b.h.p + s.h.p * y
    ((1/(2*pi)) * sqrt(1/(1 - (rho[2]^2))) 
             * exp( -0.5*((g^2) - (2*g*rho[2]*(h.p-gamma.p)) + ((h.p-gamma.p)^2))/(1-rho[2]^2) )) * s.h.p
  }
  
  GL.modified.1 <- GL_modified_nodes_weights(delt, d, 
           no.nodes.GL)
  GL.modified.2 <- GL_modified_nodes_weights(-1, 1, no.nodes.GL)
  GL.modified.3 <- GL_modified_nodes_weights(-d, -delt, 
           no.nodes.GL)
  
  term1 <- 0  
  for(j in 1:no.nodes.GL){
    hp.j <- GL.modified.1$GL.modified.nodes[j]
    w.hp.j <- GL.modified.1$GL.modified.weights[j]
    
    for (i in 1:no.nodes.GL) {
      y.i <- GL.modified.2$GL.modified.nodes[i]
      w.y.i <- GL.modified.2$GL.modified.weights[i]
      term1.val <- fn(hp.j,  y.i, rho, gamma.p, alpha, d, t.vec, 
               b.knots.vec, s.knots.vec)
      term1 <- term1 + w.y.i * w.hp.j * term1.val
    }
  }
  
  term2 <- 0
  for(j in 1:no.nodes.GL){
    hp.j <- GL.modified.3$GL.modified.nodes[j]
    w.hp.j <- GL.modified.3$GL.modified.weights[j]
    
    for (i in 1:no.nodes.GL) {
      y.i <- GL.modified.2$GL.modified.nodes[i]
      w.y.i <- GL.modified.2$GL.modified.weights[i]
      term2.val <- fn(hp.j,  y.i, rho, gamma.p, alpha, d, t.vec, 
               b.knots.vec, s.knots.vec)
      term2 <- term2 + w.y.i * w.hp.j * term2.val
    }
  }
  term1 + term2
  
}