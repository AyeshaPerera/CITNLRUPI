b_fn <- function(t, b.knots.vec, t.vec, d){
  # This module evaluates the function b(t) for
  # given (b1, b2, ..., bk). 
  #
  # b(t) =
  #  b1 * (sinc(t1 - 1) - sinc(t1 + 1))+...+bk *
  #          (sinc(tk - k) - sinc(tk + k))
  #
  # Inputs
  # t: value of the test statistic
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # t.vec: the vector (t0, t1, ..., tk)
  # d: a positive number that specifies 
  #    the interval [-d,d]
  #
  # Output
  # Value of b(t)
  #
  # Written by A.Perera 
  
  if (abs(t) >= d){
    return(0)
  }
  
  b.t.vec <- rep(0, length(b.knots.vec))
  
  for(j in 2:length(t.vec)){
    i <- t.vec[j]
    b.t.vec[j-1] <- sinc_fn((t - i) / (d / length(t.vec)) ) - sinc_fn((t + i) / (d / length(t.vec)) )
  }
  
  sum(b.knots.vec * b.t.vec) 
}

