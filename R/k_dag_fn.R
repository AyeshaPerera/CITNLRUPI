k_dag_fn <- function(h.p.1, h.p, rho, d, alpha, gamma.vec){
  # This module computes
  # k+(h, gamma, rho) = Psi(-z(1-alpha/2), z(1-alpha/2); 
  #         rho(h - gamma), 1 - rho^2)
  # 
  # Inputs
  # h.p.1: a value (h{p-1} of h vector)
  # h.p: a value (h{p} of h vector)
  # rho: vector of (rho{p-2}, rho{p})
  # d: a positive number that specifies 
  #    the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  #
  # Output
  # Value of k+(h, gamma, rho)
  #
  # Written by A.Perera May 2023
  
  h <- c(h.p.1, h.p)
  z.val <- qnorm(1 - alpha/2)
  rho.col <- matrix(rho)
  h.col <- matrix(h)
  gamma.vec.col <- matrix(gamma.vec)
  
  mu <- t(rho.col) %*% (h.col - gamma.vec.col)
  variance <- 1 - t(rho.col) %*% rho.col
  
  Psi(-z.val, z.val, mu, variance)
}
