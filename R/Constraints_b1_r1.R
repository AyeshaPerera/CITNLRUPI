Constraints_b1_r1 <- function(b1.r1.vec, t.vec, d, alpha, rho, gamma.p.vec, gamma.pm1.vec, no.nodes.GL, b.knots.vec, 
         s.knots.vec, delta, u1, u){  
  # Constraint function for the slsqp function 
  #
  # Input:
  # b1.r1.vec: the vector (b1, b2, ..., bk, r10, r11, ..., r1k)
  # t.vec: the vector (t0, t1, ..., tk)
  # d: a positive number that specifies the interval [-d,d]
  # alpha: The desired minimum coverage of the CI is 1 - alpha
  # rho: the vector of (rho_{p-1}, rho_{p}) 
  # gamma.p.vec: a vector gamma_{p} values
  # gamma.pm1.vec: a vector of gamma_{p-1} values
  # no.nodes.GL: Number of nodes used in Gauss Legendre
  # s.knots.vec: values of sc(x) at the knots
  # b.knots.vec: values of b(x) at the knots 
  # delta: a value
  # u1: a value where u1 >= SEL(gamma.pm1,0)
  # u: a value where u1 >= SEL(gamma.pm1, gamma.p)
  #
  # output:
  # Constraints
  #
  # Written by A.Perera 2023
  
  b1.r1.split <- (length(b1.r1.vec) + 1) / 2
  r.knots.vec <- b1.r1.vec[(b1.r1.split):length(b1.r1.vec)]
  b1.knots.vec <- b1.r1.vec[1:(b1.r1.split - 1)]  
  cp.vec <- CP_T2_GL(rho, d, alpha, gamma.pm1.vec, gamma.p.vec, 
           t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, 
           r.knots.vec, delta, no.nodes.GL)
  
  sel.mat <- matrix(nrow = length(gamma.pm1.vec),
           ncol = length(gamma.p.vec))
  for (i in 1:length(gamma.pm1.vec)) {
    for (j in 1:length(gamma.p.vec)) {
      gamma.vec <- c(gamma.pm1.vec[i], gamma.p.vec[j])
      sel.mat[i,j] <- SEL_GL(r.knots.vec, s.knots.vec, t.vec, 
               alpha, d, gamma.vec, no.nodes.GL, delta)
    }
  }  
  pos.gamma.pm1.vec <- gamma.pm1.vec[gamma.pm1.vec>=0]
  sel.vec <- rep(0, length(pos.gamma.pm1.vec))
  for (j in 1:length(pos.gamma.pm1.vec)) {
    gamma.vec <- c(pos.gamma.pm1.vec[j], 0)
    sel.vec[j] <- SEL_GL(r.knots.vec, s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta) 
  }
  step_len <- d/length(t.vec)
  t.values <- seq(0, (d-step_len), by = step_len)  
  r1.vec <- rep(0, length(t.values))
  for(t in 1:length(t.values)) {
  	r1.vec <- r1_fn(t, r.knots.vec, t.vec, alpha, d)
  }
  epsilon = -qnorm(1-alpha/2)/0.5
  cp.minus.vec <- cp.vec - (1-alpha)
  sel.gammapm1.0 <- u1 - sel.vec
  max.sel.mat <- u - max(sel.mat)
  h.vec <- c(sel.gammapm1.0, cp.minus.vec, max.sel.mat, 
           (r1.vec - epsilon))  
  return(h.vec)  
}