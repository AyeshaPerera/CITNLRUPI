slsqp_option <- function (x0, fn, gr = NULL, lower = NULL, upper = NULL, hin = NULL,
                          hinjac = NULL, heq = NULL, heqjac = NULL, nl.info = FALSE,
                          control = list(), ...)
{
  opts <- nl.opts(control)
  opts["algorithm"] <- "NLOPT_LD_SLSQP"
  opts["print_level"] <- 3
  fun <- match.fun(fn)
  fn <- function(x) fun(x, ...)
  if (is.null(gr)) {
    gr <- function(x) nl.grad(x, fn)
  }
  else {
    .gr <- match.fun(gr)
    gr <- function(x) .gr(x, ...)
  }
  if (!is.null(hin)) {
    if (getOption("nloptr.show.inequality.warning")) {
      message("For consistency with the rest of the package the inequality sign may be switched from >= to <= in a future nloptr version.")
    }
    .hin <- match.fun(hin)
    hin <- function(x) (-1) * .hin(x)
    if (is.null(hinjac)) {
      hinjac <- function(x) nl.jacobian(x, hin)
    }
    else {
      .hinjac <- match.fun(hinjac)
      hinjac <- function(x) (-1) * .hinjac(x)
    }
  }
  if (!is.null(heq)) {
    .heq <- match.fun(heq)
    heq <- function(x) .heq(x)
    if (is.null(heqjac)) {
      heqjac <- function(x) nl.jacobian(x, heq)
    }
    else {
      .heqjac <- match.fun(heqjac)
      heqjac <- function(x) .heqjac(x)
    }
  }
  s.time <- Sys.time()
  S0 <- nloptr(x0, eval_f = fn, eval_grad_f = gr, lb = lower,
               ub = upper, eval_g_ineq = hin, eval_jac_g_ineq = hinjac,
               eval_g_eq = heq, eval_jac_g_eq = heqjac, opts = opts)
  e.time <- Sys.time()
  SO.out <- list()
  SO.out[["SO"]] = S0
  SO.out[["time"]] = e.time - s.time
  if (nl.info)
    print(SO.out)
  S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
             convergence = S0$status, message = S0$message)
  return(S1)
}