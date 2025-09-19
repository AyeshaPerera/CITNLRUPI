s_fn <- function(t, s.knots.vec, t.vec, alpha, d){
  # This module evaluates the function s(t) for
  # given (s0, s1, ..., sk). 
  #
  # s(t) = qnorm(1-alpha/2) + s_c(t)
  #
  # Where s_c(t) = 
  #  s0 * sinc(t0) + s1 * (sinc(t1 - 1) + sinc(t1 + 1)) + ... + 
  #  sk * (sinc(tk - k) + sinc(tk + k))
  #
  # Inputs :
  # t: value of the test statistic
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  #
  # Output
  # Value of s(t)
  #
  # Written by A.Perera
  
  z <- qnorm((1 - alpha/2))
  
  if(abs(t) >= d){
    return(z)
  }
  
  len.t.vec <- length(t.vec)
  
  s.t.vec <- rep(0, len.t.vec)
  s.t.vec[1] <- sinc_fn(t / (d / len.t.vec))
  
  for (j in 2:len.t.vec) {
    i <- t.vec[j]
    s.t.vec[j] <- (sinc_fn( (t + i) / (d / length(t.vec)) )+ sinc_fn((t - i) / (d / length(t.vec)) ) )
  }
  
  s.c <- sum(s.knots.vec * s.t.vec)
  s.t <- z + s.c
  s.t
}


