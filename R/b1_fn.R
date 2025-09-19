b1_fn <- function(x, b1.knots.vec, t.vec, d){
  # This module evaluates the function b1(x) for
  # given (b1, b2, ..., bk). 
  #
  # b1(x) = b11 * (sinc(x1 - 1) - sinc(x1 + 1)) + ... + 
  #         b1k * (sinc(xk - k) - sinc(xk + k))
  #
  # Inputs
  # t: value of the test statistic
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # t.vec: the vector (t0, t1, ..., tk)
  # d: a positive number that specifies 
  #    the interval [-d,d]
  #
  # Output
  # Value of b1(x)
  #
  # Written by A.Perera 
  
  if (abs(x) >= d){
    return(0)
  }
  
  b1.x.vec <- rep(0, length(b1.knots.vec))
  
  for(j in 2:length(t.vec)){
    i <- t.vec[j]
    b1.x.vec[j-1] <- sinc_fn((x - i) / (d / length(t.vec)) ) - sinc_fn((x + i) / (d / length(t.vec)) ) 
  }
  
  sum(b1.knots.vec * b1.x.vec) 
}