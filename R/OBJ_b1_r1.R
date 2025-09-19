OBJ_b1_r1 <- function(b1.r1.vec, gamma.pm1.vec, s.knots.vec, t.vec, alpha, d, no.nodes.GL, delta){
  # Computes the values of 
  # l(max_(gamma_{p-1} >= 0) SEL(gamma_{p-1},0) 
  #         + (1 - l) SEL(0,0)
  #
  # Input:
  # b1.r1.vec: the vector (b1, b2, ..., bk, r10, r11, ..., r1k)
  # gamma.pm1.vec: a vector of gamma_{p-1} values
  # s.knots.vec: values of sc(x) at the knots
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The desired minimum coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # no.nodes.GL: number of nodes used in Gauss Legendre
  # delta: a value
  #
  # output:
  # value of l(max_(gamma_{p-1} >= 0) SEL(gamma_{p-1},0) 
  #         + (1 - l) SEL(0,0)
  #
  # Written by A.Perera 2023 
  l <- 0.03
  b1.r1.split <- (length(b1.r1.vec) + 1) / 2
  r.knots.vec <- b1.r1.vec[(b1.r1.split):length(b1.r1.vec)]
  
  sel.0.0 <- SEL_GL(r.knots.vec, s.knots.vec, t.vec, alpha, d, 
           c(0,0), no.nodes.GL, delta)
  pos.gamma.pm1.vec <- gamma.pm1.vec[gamma.pm1.vec>=0]
  sel.vec <- rep(0, length(pos.gamma.pm1.vec))
  for (j in 1:length(pos.gamma.pm1.vec)) {
    gamma.vec <- c(pos.gamma.pm1.vec[j], 0)
    sel.vec[j] <- SEL_GL(r.knots.vec, s.knots.vec, t.vec, alpha, 
             d, gamma.vec, no.nodes.GL, delta)
  }
  l*max(sel.vec) + (1-l) * sel.0.0
}
