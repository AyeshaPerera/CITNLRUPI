r1_fn <- function(t, r.knots.vec, t.vec, alpha, d){
  # This module evaluates the function r1(t) for
  # given (r10, r11, ..., r1k). 
  #
  # r1(t) = r10 * sinc(t0) + r11 * (sinc(t1 - 1) + sinc(t1 + 1)) 
  #         + ... + r1k * (sinc(tk - k) + sinc(tk + k))
  #
  # Inputs :
  # t: value of the test statistic
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  #
  # Output
  # Value of s(t)
  #
  # Written by A.Perera 
  
  len.t.vec <- length(t.vec)
  
  if (abs(t) >= d){
    return(0)
  }
  
  r1.t.vec <- rep(0, len.t.vec)
  r1.t.vec[1] <- sinc_fn(t / (d / len.t.vec))
  
  for (j in 2:len.t.vec) {
    i <- t.vec[j]
    r1.t.vec[j] <- (sinc_fn( (t + i) / (d / length(t.vec)) )+ sinc_fn((t - i) / (d / length(t.vec)) ) )
  }
  
  sum(r.knots.vec * r1.t.vec)
  
}


