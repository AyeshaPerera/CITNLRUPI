CP_T2_GL_t3 <- function(rho, gam, alpha, d, t.vec, b.knots.vec, s.knots.vec, delt, no.nodes.GL){
  # This module computes
  #    delta       1
  #  integral integral  (f_{G, H_p}(b(h_p) + s(h_p)y, h_p) 
  #     -delta       -1
  #			* s_1(b(h_p) + s(h_p)y, h_p) * s(h_p) dy dh_p)
  # using the Gauss Legendre quadrature function
  #
  # Input
  # rho: vector of (rho{p-2}, rho{p})
  # gam: vector of (gamma{p-1}, gamma{p})
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # delt: tuning constant for w1(x) were 0 < delt <= d 
  #
  # Output
  #    delta       1
  #  integral integral  (f_{G, H_p}(b(h_p) + s(h_p)y, h_p) 
  #     -delta       -1
  #			* s_1(b(h_p) + s(h_p)y, h_p) * s(h_p) dy dh_p)
  # Written by A.Perera June 2023  
  trans_fn <- function(h.p, y, rho, gam, alpha, d, t.vec, 
           b.knots.vec, s.knots.vec){
    b.h.p <- b_fn(t=h.p, b.knots.vec=b.knots.vec, t.vec=t.vec, 
             d=d)
    s.h.p <- s_fn(t=h.p, s.knots.vec=s.knots.vec, t.vec=t.vec, 
             alpha = alpha, d=d)
    g <- b.h.p + s.h.p * y
    ((1/(2*pi)) * sqrt(1/(1 - (rho[2]^2))) 
       * exp( -0.5*((g^2) - (2*g*rho[2]*(h.p-gam[2])) + ((h.p-gam[2])^2))/(1-rho[2]^2) )) * s1_fn(h.p, g, rho, gam, d) * s.h.p
  }  
  GL.modified.1 <- GL_modified_nodes_weights(-delt, delt,
           no.nodes.GL)
  GL.modified.2 <- GL_modified_nodes_weights(-1, 1, no.nodes.GL)  
  term <- 0
  for(j in 1:no.nodes.GL){
    hp.j <- GL.modified.1$GL.modified.nodes[j]
    w.hp.j <- GL.modified.1$GL.modified.weights[j]
    for (i in 1:no.nodes.GL) {
      y.i <- GL.modified.2$GL.modified.nodes[i]
      w.y.i <- GL.modified.2$GL.modified.weights[i]
      term.val <- trans_fn(hp.j,  y.i, rho, gam, alpha, d, 
               t.vec, b.knots.vec, s.knots.vec)
      term <- term + w.y.i * w.hp.j * term.val
    }
  }
  term
}
