SEL_GL_t1 <- function(s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL){
  # This module computes
  #     1           d
  # ------------ integrate (s(h{p})) * psi(h{p} - gamma{p}) dh{p}
  # z(1-alpha/2)   -d
  #
  # Inputs
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # no.nodes.GL: the number of Gauss Legendre nodes
  #
  # Output
  # first term integration part of the equation to scaled expected length
  #
  # Written by A.Perera May 2023
  
  z <- qnorm(1-alpha/2)
  
  fn <- function(h.p, s.knots.vec, t.vec, alpha, d, gamma.p){
    # This function represents (s(h{p})) * psi(h{p} - gamma{p})
    s.h.p <- s_fn(h.p, s.knots.vec, t.vec, alpha, d)
    (s.h.p) * dnorm(h.p - gamma.p)
  }
  
  GL.modified <- GL_modified_nodes_weights(-d, d, no.nodes.GL)
  sel <- 0
  for (i in 1:no.nodes.GL) {
    h.p.i <- GL.modified$GL.modified.nodes[i]
    w.h.p.i <- GL.modified$GL.modified.weights[i]
    sel.val <- fn(h.p.i, s.knots.vec, t.vec, alpha, d, gamma.vec[2])
    sel <- sel + w.h.p.i * sel.val 
  }
  sel / z
}