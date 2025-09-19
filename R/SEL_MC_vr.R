#' Compute the scaled expected length of the CIUUPI for a given vector of
#' \eqn{\gamma} using Monte Carlo simulation with variance reduction
#'
#' This function computes the scaled expected length using Monte Carlo
#' simulation with variance reduction.
#' The scaled expected length is
#' \deqn{
#' {\rm SEL}_{\rm MC}(\bm{\gamma}) = {\rm E}\left(\frac{s(H_p) + w_1(|H_p|) r_1(H_{p-1})}
#' { z_{1-\alpha/2}}\right).
#' }
#' Here \eqn{\bm{H} \sim N(\bm{\gamma}, \bm{I})} where
#' \eqn{\bm{H} = (H_{p-1}, H_p)}.
#' We estimate \eqn{{\rm SEL}_{\rm MC}(\bm{\gamma})} by
#' \deqn{
#' \frac{1}{N} \sum_{r=1}^{N} \left(\frac{s(h_p^r) + w_1(|h_p^r|)
#' r_1(h_{p-1}^r)}{ z_{1-\alpha/2}}\right),
#' }
#' where \eqn{N} is the number of independent simulation runs and
#' \eqn{(h_{p-1}^r, h_p^r)} is the observation of \eqn{\bm{H}} in the
#' \eqn{r^{th}} simulation run.
#'
#' @param gamma.vec Vector of \eqn{(\gamma_{p-1}, \gamma_{p})}
#' @param N Number of independent simulation runs
#' @param s.knots.vec The vector \eqn{(s_0, s_1, \dots, s_k)}
#' @param r.knots.vec The vector \eqn{(r_0, r_1, \dots, r_k)}
#' @param t.vec The vector \eqn{(t_0, t_1, \dots, t_k)}
#' @param alpha The nominal coverage of the CI is \eqn{1 - \alpha}
#' @param d A positive number that specifies the interval \eqn{[-d,d]}
#' @param delta A positive number that specifies the interval \eqn{[0,d]}
#'
#' @returns
#' The scaled expected length computed using Monte Carlo simulation with
#' variance reduction.
#' @export
#'
SEL_MC_vr <- function(gamma.vec, N, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta){

  set.seed(100)

  # This module calculates the scaled expected length using Monte Carlo simulation with common random numbers
  #
  # Input
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # N: Number of independent simulation runs
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # delta: tuning constant for w1(x) were 0 < delta <= d
  #
  # Output:
  # Scaled expected length estimate
  #
  # Written by A.Perera June 2023

  # List to store the estimated scaled expected length and it's standard error
  sel.list <- list()

  # generate random values
  values <- rmvnorm(N*length(gamma.vec), mean = gamma.vec, sigma = matrix(c(1, 0, 0, 1), nrow = 2))

  func <- function(h, s.knots.vec, r.knots.vec, t.vec, alpha, d){

    z <- qnorm(1 - alpha/2)
    h.p <- h[2]
    h.p.1 <- h[1]
    s.h.p <- s_fn(h.p, s.knots.vec, t.vec, alpha, d)
    w.1.h.p <- w1_fn(abs(h.p), d, delta)
    r.1.h.p <- r1_fn(h.p.1, r.knots.vec, t.vec, alpha, d)
    (s.h.p + w.1.h.p * r.1.h.p) / z
  }

  estimates <- apply(values, 1, func, s.knots.vec=s.knots.vec, r.knots.vec=r.knots.vec, t.vec=t.vec, alpha=alpha, d=d)
  mc.estimate <- mean(estimates)
  s <- sqrt(sum((estimates - mc.estimate)^2) / (N - 1))
  se <- s / sqrt(N)
  sel.list$Estimate <- mc.estimate
  sel.list$SE <- se
  sel.list
}
