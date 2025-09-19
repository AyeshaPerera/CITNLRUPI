s1_fn<- function(h.p, g, rho, gamma.vec, d){
  # This module computes s1(g, h_p)
  # = integral    f_{H_{p-1}|G, H_p}(h_{p-1}|g, h_p) dh_{p-1}
  #   |h_{p-1}|>d 
  # using the integrate() function
  #
  # Input
  # h.p: a value
  # g: a value
  # rho: vector of (rho{p-1}, rho{p})
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # t.knots.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # delta: tuning constant for w1(x) were 0 < delta <= d 
  #
  # Output
  #    s_1(g, h_p)
  #
  # Written by A.Perera June 2023
  
  mu <- gamma.vec[1] + (1 / (1 - rho[2]^2)) * 
           (g * rho[1] - rho[2] * rho[1] * (h.p - gamma.vec[2]))
  Sigma <- (1 / (1 - rho[2]^2)) * (1 - rho[1]^2 - rho[2]^2)
  pnorm(-d, mean = mu, sd = sqrt(Sigma)) + 1 -
           pnorm(d, mean = mu, sd = sqrt(Sigma))

}
