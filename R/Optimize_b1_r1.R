Optimize_b1_r1 <- function(start.vec, gamma.pm1.vec, 
         gamma.p.vec, rho, b.knots.vec, s.knots.vec, t.vec, alpha, d, no.nodes.GL, delta, u1, u){
  # Optimization
  #
  # Input:
  # start.vec: the start vector
  # gamma.pm1.vec: a vector of gamma(p-1) values
  # gamma.p.vec: a vector of gamma(p) values
  # rho: the vector of (rho_{p-1}, rho_{p}) 
  # b.knots.vec: values of b(x) at the knots
  # s.knots.vec: values of sc(x) at the knots
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The desired minimum coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # no.nodes.GL: Number of nodes used in Gauss Legendre
  # delta: a value
  # u1: a value where u1 >= sel(gamma.pm1,0)
  # u: a value where u1 >= SEL(gamma.pm1, gamma.p)
  #
  # output:
  # Constraints
  #
  # Written by A.Perera 2023
  
  
  obj_fn <- functional::Curry(OBJ_b1_r1,
           gamma.pm1.vec=gamma.pm1.vec, s.knots.vec=s.knots.vec, 
           t.vec=t.vec, alpha=alpha, d=d, 
           no.nodes.GL=no.nodes.GL, delta=delta)
  
  constraints_fn <- functional::Curry(Constraints_b1_r1,
           t.vec=t.vec, d=d, alpha=alpha, rho=rho, 
           gamma.p.vec=gamma.p.vec, gamma.pm1.vec=gamma.pm1.vec, 
           no.nodes.GL=no.nodes.GL, b.knots.vec=b.knots.vec, 
           s.knots.vec=s.knots.vec, delta=delta, u1=u1, u=u)
  
  up <- rep(2, length(start.vec))
  low <- c(rep(-2, 5), rep(-0.5 * qnorm(1 - alpha/2), 6))
  
  sl <- slsqp_option(x0=start.vec, fn=obj_fn , 
           hin = constraints_fn, 
           control = list(xtol_rel = 1e-6), nl.info = TRUE, 
           lower = low, upper = up)
  
  sl
}