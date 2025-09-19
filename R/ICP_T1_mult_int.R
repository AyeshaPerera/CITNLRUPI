ICP_T1_mult_int <- function(rho, d, alpha, gamma.vec, t.vec, 
         b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, 
         delta, no.nodes.GL){
  # This module computes 
  #    d        Inf
  #  integral integral  (k(h, gamma, rho) - k+(h, gamma, rho))
  #     -d       -Inf
  #         * phi(h - gamma) dh{p-1} dh{p}
  # using multiple integrate functions
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

  integral1.1 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = d, 
                      upper = 2*d)$value
  }
  term.2.1.1 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.1.1 <- term.2.1.1 + w.h.p.j * integral1.1(h.p.j)
    
  }
  
  integral1.2 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = 2*d,  
                      upper = 3*d)$value
  }
  term.2.1.2 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.1.2 <- term.2.1.2 + w.h.p.j * integral1.2(h.p.j)
    
  }
  
  integral1.3 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = 3*d, 
                      upper = 4*d)$value
  }
  term.2.1.3 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.1.3 <- term.2.1.3 + w.h.p.j * integral1.3(h.p.j)
    
  }
  
  integral1.4 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = 4*d, 
                      upper = 5*d)$value
  }
  
  term.2.1.4 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.1.4 <- term.2.1.4 + w.h.p.j * integral1.4(h.p.j)
    
  }
  
  integral1.5 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = 5*d, 
                      upper = Inf)$value
  }
  term.2.1.5 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.1.5 <- term.2.1.5 + w.h.p.j * integral1.5(h.p.j)
    
  }
  
  term.2.1 <- term.2.1.1 + term.2.1.2 + term.2.1.3 + term.2.1.4 
           + term.2.1.5
  
  integral2.1 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = -2*d, 
                      upper = -d)$value
  }
  term.2.2.1 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.2.1 <- term.2.2.1 + w.h.p.j * integral2.1(h.p.j)
    
  }
  
  integral2.2 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = -3*d, 
                      upper = -2*d)$value
  }
  term.2.2.2 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.2.2 <- term.2.2.2 + w.h.p.j * integral2.2(h.p.j)
    
  }
  
  integral2.3 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = -4*d, 
                      upper = -3*d)$value
  }
  term.2.2.3 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.2.3 <- term.2.2.3 + w.h.p.j * integral2.3(h.p.j)
    
  }
  
  integral2.4 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = -5*d, 
                      upper = -4*d)$value
  }
  term.2.2.4 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.2.4 <- term.2.2.4 + w.h.p.j * integral2.4(h.p.j)
    
  }
  
  integral2.5 <- function(h.p) {
    integrate(function(h.p.1){
             sapply(h.p.1, I_fn,h.p=h.p, rho=rho, d=d, 
                      alpha=alpha, gamma.vec=gamma.vec, 
                      t.vec=t.vec, delta=delta, 
                      b.knots.vec=b.knots.vec, 
                      b1.knots.vec=b1.knots.vec, 
                      s.knots.vec=s.knots.vec, 
                      r.knots.vec=r.knots.vec)},lower = -Inf, 
                      upper = -5*d)$value
  }
  term.2.2.5 <- 0
  for(j in 1:no.nodes.GL){
    h.p.j <- GL.modified$GL.modified.nodes[j]
    w.h.p.j <- GL.modified$GL.modified.weights[j]
    term.2.2.5 <- term.2.2.5 + w.h.p.j * integral2.5(h.p.j)
    
  }
  
  term.2.2 <- term.2.2.1 + term.2.2.2 + term.2.2.3 + term.2.2.4 + term.2.2.5
  term.1 + term.2.1 + term.2.2
}
