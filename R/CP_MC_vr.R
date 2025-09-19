#' Compute the coverage probability of the CIUUPI for a given vector of
#' \eqn{\gamma} and \eqn{\rho} using Monte Carlo simulation with variance reduction
#'
#' This function computes the coverage probability using Monte Carlo
#' simulation with variance reduction.
#' The coverage probability is
#' \deqn{
#' {\rm CP}_{\rm MC}(\bm{\gamma}, \bm{\rho})
#' = \int k(\bm{h}, \bm{\gamma}, \bm{\rho}) \,\bm{\phi}(\bm{h}-\bm{\gamma}) \, d\bm{h}
#' \\[6pt]
#' = \text{E}\!\left[ k(\bm{H}, \bm{\gamma}, \bm{\rho}) \right]
#' }
#' Here \eqn{\bm{H} \sim N(\bm{\gamma}, \bm{I})} where
#' \eqn{\bm{H} = (H_{p-1}, H_p)}.
#' We estimate \eqn{{\rm CP}_{\rm MC}(\bm{\gamma}, \bm{\rho})} by
#' \deqn{
#' \frac{1}{N} \sum_{r=1}^{N} k(\bm{h}^r, \bm{\gamma}, \bm{\rho}),
#' }
#' where \eqn{N} is the number of independent simulation runs and
#' \eqn{(h_{p-1}^r, h_p^r)} is the observation of \eqn{\bm{H}} in the
#' \eqn{r^{th}} simulation run.
#'
#' @param gamma.vec Vector of \eqn{(\gamma_{p-1}, \gamma_{p})}
#' @param rho Vector of \eqn{(\rho_{p-1}, \rho_{p})}
#' @param N Number of independent simulation runs
#' @param b.knots.vec The vector \eqn{(b(1), b(2), \dots, b(k))}
#' @param b1.knots.vec The vector \eqn{(b_1(0), b_1(1), \dots, b_1(k))}
#' @param s.knots.vec The vector \eqn{(s(0), s(1), \dots, s(k))}
#' @param r.knots.vec The vector \eqn{(r_1(0), r_1(1), \dots, r_1(k))}
#' @param t.vec The vector \eqn{(t_0, t_1, \dots, t_k)}
#' @param alpha The nominal coverage of the CI is \eqn{1 - \alpha}
#' @param d A positive number that specifies the interval \eqn{[-d,d]}
#' @param delta A positive number that specifies the interval \eqn{[0,d]}
#'
#' @returns
#' The coverage probability computed using Monte Carlo simulation with
#' variance reduction.
#' @export
#'
#'
CP_MC_vr <- function(gamma.vec, rho, N, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta){

  set.seed(890)
  # This module calculates the coverage probability using Monte
  # Carlo simulation variance reduction by conditioning
  #
  # Input
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # rho: vector of (rho{p-1}, rho{p})
  # N: Number of independent simulation runs
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies the interval [-d,d]
  # delta: tuning constant for w1(x) were 0 < delta <= d
  #
  # Output:
  # coverage probability estimate
  #
  # Written by A.Perera June 2023

  # Simulate random vector [H] using multivariate normal distribution which (gamma, I)
  obs.h <- rmvnorm(N, mean = gamma.vec, sigma = matrix(c(1, 0, 0, 1), nrow = 2))

  k.obs.h <- rep(0, N)
  for (i in 1:N) {
    k.obs.h[i] <- k_fn(obs.h[i, 1], obs.h[i, 2], rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta)
  }

  cp.estimate <- mean(k.obs.h)
  s.sq <- sum((k.obs.h - cp.estimate)^2) / (N - 1)
  se <- sqrt(s.sq / N)

  # List to store the estimated coverage probability and it's standard error
  list(Estimate=cp.estimate, SE=se)

}
