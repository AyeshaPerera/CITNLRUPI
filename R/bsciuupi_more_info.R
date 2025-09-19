library(ciuupi)


bsciuupi_more_info <- function (alpha, rho , natural = 1, a=NULL, c=NULL, x=NULL){
  gams <- seq(0, 8, by = 0.05)
  n.iter <- 5
  d <- 6
  n.ints <- 6
  n.nodes <- 5
  if (is.null(rho)) {
    qrstr <- qr(x)
    R <- qr.R(qrstr)
    XTXinv <- chol2inv(R)
    rho <- (t(a) %*% XTXinv %*% c)/sqrt(t(a) %*% XTXinv %*% 
                                          a %*% t(c) %*% XTXinv %*% c)
    rho <- as.numeric(rho)
  }
  cat("Computing the vector (b(1),...,b(5),s(0),...,s(5)) that specifies the CIUUPI... ")
  if (rho == 0) {
    new.par <- ciuupi:::standard_CI(d, n.ints, alpha)
  }else {
    lambda <- ciuupi:::compute_lambda(rho, alpha, n.iter, d, n.ints, 
                                      n.nodes, gams, natural)
    new.par <- optimize_knots_info(lambda, rho, alpha, gams, 
                                   d, n.ints, n.nodes, natural)
  }
  cat("DONE", "\n")
  out <- list("new.par"=new.par, "lambda"=lambda)
  out
}

