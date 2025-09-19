Psi <- function (x, y, mu , variance){
  # This function calculates
  # Psi(x, y, mu, variance) which equals to
  # P(x le Z le y) = P(Z le y) - P(Z le x)
  # where Z ~ N(mu, variance).
  #
  # Inputs
  # x: given value
  # y: given value that is greater than or equal
  #    to x
  # mu: mean of the normal distribution
  # variance : variance of the normal distribution
  #
  # Output
  # Psi(x, y, mu, variance)
  #
  # Written by N. Ranathunga in September 2020
  
  sigma <- sqrt(variance)
  term1 <- pnorm(y, mean = mu , sd = sigma)
  term2 <- pnorm(x, mean = mu , sd = sigma)
  
  term1 - term2
  
}