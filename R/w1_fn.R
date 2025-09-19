w1_fn <- function(x, d, delta){
  # This module computes 
  # w1(x) = 
  #  exp(2) * exp(- 2 / (1 - (x / delta)^2)) for x < delta
  #  0 for x >= d
  #
  # Inputs
  # x: a number
  # d: a tuning constant
  # delta: tuning constant for w1(x) were 0 < delta <= d. 
  #
  # Output
  # Value of w1(x)
  #
  # Written by A.Perera May 2023
  
  if(abs(x) < d && abs(x) < delta){
    out <- exp(2) * exp(-2 / (1 - (x / delta)^2))
  }else{
    out <- 0
  }
  out
  
}