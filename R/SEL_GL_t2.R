SEL_GL_t2 <- function(r.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta){
  # This module computes
  #     1           d
  # ------------ (integrate (w1(|h{p}|))psi(h{p} - gamma{p}) dh{p} 
  # z(1-alpha/2)   -d
  #     d
  # integrate r1(h{p-1}) psi(h{p-1} - gamma{p-1}) dh{p-1})
  #   -d
  #
  # Inputs
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # no.nodes.GL: the number of Gauss Legendre nodes
  # delta: tuning constant for w1(x) were 0 < delta <= d 
  #
  # Output
  # second term integration part of the equation to scaled expected length
  #
  # Written by A.Perera May 2023
  
  z <- qnorm(1-alpha/2)
  
  fn1 <- function(h.p, d, gamma.p, delta){
    # This function represents (w1(|h{p}|))psi(h{p} - gamma{p})
    w1_fn(abs(h.p), d, delta) * dnorm(h.p - gamma.p)
  }
  
  GL.modified.1 <- GL_modified_nodes_weights(-delta, delta, no.nodes.GL)
  sel1 <- 0
  for (i in 1:no.nodes.GL) {
    h.p.i <- GL.modified.1$GL.modified.nodes[i]
    w.h.p.i <- GL.modified.1$GL.modified.weights[i]
    sel1.val <- fn1(h.p.i, d, gamma.vec[2], delta)
    sel1 <- sel1 + w.h.p.i * sel1.val 
  }
  
   
  fn2 <- function(h.p.1, r.knots.vec, t.vec, alpha, d, gamma.p.1){
    # This function represents r1(h{p-1}) psi(h{p-1} - gamma{p-1})
    r1_fn(h.p.1, r.knots.vec, t.vec, alpha, d) * dnorm(h.p.1 - gamma.p.1)
  }
  
  
  GL.modified.2 <- GL_modified_nodes_weights(-d, d, no.nodes.GL)
  sel2 <- 0
  for (j in 1:no.nodes.GL) {
    y.j <- GL.modified.2$GL.modified.nodes[j]
    w.y.j <- GL.modified.2$GL.modified.weights[j]
    sel2.val <- fn2(y.j, r.knots.vec, t.vec, alpha, d, gamma.vec[1])
    sel2 <- sel2 + w.y.j * sel2.val 
  }
  sel1 * sel2/ z
}
