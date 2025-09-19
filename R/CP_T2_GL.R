CP_T2_GL <- function(rho, d, alpha, gamma.pm1.vec, gamma.p.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta, no.nodes.GL){
  
  # Computes coverage probability using Gauss Legendre quadrature
  #
  # Input
  # rho: vector of (rho{p-1}, rho{p})
  # d: a positive number that specifies the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.pm1.vec: vector of gamma{p-1}
  # gamma.p.vec: vector of gamma{p}
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # delta: tuning constant for w1(x) were 0 < delta <= d 
  # no.nodes.GL: Number of nodes
  #
  # Output
  # Coverage probability
  #
  # Written by A. Perera June 2023
  
  z <- qnorm(1 - alpha / 2)
  Sigma <- matrix(c(1, rho[2], rho[2], 1), nrow=2)
  
  term1.4.vec <- rep(0, length(gamma.p.vec))
  for (k in 1:length(gamma.p.vec)) {
    mu <- c(0, gamma.p.vec[k])
    term1 <- pmvnorm(lower=c(-z, -d), upper=c(z, d), algorithm = "Miwa", mean=mu, sigma=Sigma)[1]
    term4 <- CP_T2_GL_t4(rho, gamma.p.vec[k], alpha, d, t.vec, b.knots.vec, s.knots.vec, delta, no.nodes.GL)
    term1.4.vec[k] <- term4 - term1
  }
  cp.mat <- matrix(nrow = length(gamma.pm1.vec), ncol = length(gamma.p.vec))

  for (i in 1:length(gamma.pm1.vec)) {
    for (j in 1:length(gamma.p.vec)) {
      gamma.vec <- c(gamma.pm1.vec[i], gamma.p.vec[j])
      term2 <- k_phi_int(rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta, no.nodes.GL)
      term3 <- CP_T2_GL_t3(rho, gamma.vec, alpha, d, t.vec, b.knots.vec, s.knots.vec, delta, no.nodes.GL)
      cp.mat[i,j] <- 1 - alpha + term1.4.vec[j] + term2 + term3
    }
  }
  cp.mat
}