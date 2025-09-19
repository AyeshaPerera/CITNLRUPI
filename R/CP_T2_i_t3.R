CP_T2_int_t3 <- function(rho, gamma.vec, alpha, d, t.vec, b.knots.vec, s.knots.vec, delta){

  # This module computes
  #    delta       1
  #  integral integral  f_{G, H_p}(b(h_p) + s(h_p)y, h_p)
  #     -delta       -1
  #         * s_1(b(h_p) + s(h_p)y, h_p) * s(h_p) dy dh_p
  # using the integrate() function
  #
  # Input
  # rho: vector of (rho{p-1}, rho{p})
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # delta: tuning constant for w1(x) were 0 < delta <= d
  #
  # Output
  #    delta       1
  #  integral integral  f_{G, H_p}(b(h_p) + s(h_p)y, h_p)
  #     -delta       -1
  #         * s_1(b(h_p) + s(h_p)y, h_p) * s(h_p) dy dh_p
  #
  # Written by A.Perera June 2023

  func <- function(y,  h.p, rho, d, alpha, gamma.vec, t.vec,
           b.knots.vec, s.knots.vec){
    b.h.p <- b_x(t=h.p, b.knots.vec=b.knots.vec, t.vec=t.vec,
             d=d)
    s.h.p <- s_x(t=h.p, s.knots.vec=s.knots.vec, t.vec=t.vec,
             alpha = alpha, d=d)
    g <- b.h.p + s.h.p * y
    X <- c(g, h.p)
    mu <- c(0, gamma.vec[2])
    Sigma <- matrix(c(1, rho[2], rho[2], 1), nrow=2)

    dmvnorm(X, mean = mu, sigma = Sigma) * s_g_hp(h.p, g, rho, gamma.vec, d) * s.h.p
  }

  integral_2 <- function(h.p) {
    integrate(function(y){sapply(y, func,h.p=h.p, rho=rho, d=d,
             alpha=alpha, gamma.vec=gamma.vec, t.vec=t.vec,
             b.knots.vec=b.knots.vec,
             s.knots.vec=s.knots.vec)},-1, 1,
             rel.tol = 10^(-10))$value
  }
  integrate(function(h.p){sapply(h.p, integral_2)}, -delta,
           delta, rel.tol = 10^(-10))$value
}
