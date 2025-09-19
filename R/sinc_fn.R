sinc_fn <- function(t){
  # This module computes
  # sinc(t) = sin(pi t) / (pi t) for t non-zero
  #         = 1 for t=0.
  # This module takes account of the difficulty
  # in accurately computing sinc(t) using this
  # formula for non-zero t very close to 0.
  #
  # Input
  # t: a number
  #
  # Output
  # sinc(t)
  #
  # Written by P.Kabaila in Dec 22
  
  term <- pi * t
  if (abs(t) < 1e-3){
    out <- 1 - (term^2)/6 + (term^4)/120
    return(out)
  }else{
    return(sin(term) / term)
  }
  
}