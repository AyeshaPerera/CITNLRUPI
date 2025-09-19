optimize_knots_info <- function (lambda, rho, alpha, gams, d, n.ints, n.nodes, natural) 
{
  c.alpha <- stats::qnorm(1 - alpha/2)
  start <- ciuupi:::standard_CI(d, n.ints, alpha)
  low <- c(rep(-100, n.ints - 1), rep(0.5, n.ints))
  up <- c(rep(100, n.ints - 1), rep(200, n.ints))
  obj_fun <- functional::Curry(ciuupi:::objective, lambda = lambda, 
                               d = d, n.ints = n.ints, n.nodes = n.nodes, alpha = alpha, 
                               natural = natural)
  cons_fun <- functional::Curry(ciuupi:::constraints_slsqp_gausslegendre, 
                                gams = gams, rho = rho, d = d, n.ints = n.ints, alpha = alpha, 
                                n.nodes = n.nodes, natural = natural)
  res <- nloptr::slsqp(start, obj_fun, hin = cons_fun, lower = low, 
                       upper = up, nl.info = TRUE)
  new.par <- res$par
  out <- new.par
}
