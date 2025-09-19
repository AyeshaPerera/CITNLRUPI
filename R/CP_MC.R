#' Compute the coverage probability of the CIUUPI for a given vector of
#' \eqn{\gamma} and \eqn{\rho} using Monte Carlo simulation
#'
#' This function computes the coverage probability using Monte Carlo
#' simulation with variance reduction.
#' The coverage probability is
#' \deqn{
#' {\rm CP}(\bm{\gamma}, \bm{\rho}) = P\big(l(\bm{H}) \leq G \leq u(\bm{H}) \big)
#' }
#' Here \eqn{\bm{H} \sim N(\bm{\gamma}, \bm{I})} where
#' \eqn{\bm{H} = (H_{p-1}, H_p)}.
#' We estimate \eqn{{\rm CP}_{\rm MC}(\bm{\gamma}, \bm{\rho})} by
#' \deqn{
#' \frac{1}{N} \, \sum_{r=1}^N \, \bm{1} \big(l(\bm{h}^r) \leq g^r \leq u(\bm{h}^r) \big),
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
#' The coverage probability computed using Monte Carlo simulation.
#' @export
#'
#'
CP_MC <- function(gamma.vec, rho, N, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, t.vec, alpha, d, delta) {

  # This module calculates the coverage probability using
  # Monte Carlo simulation
  #
  # Inputs
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

  # Define parameters for multivariate normal distribution
  mean.vec <- c(0, gamma.vec)
  cov.mat <- matrix(c(1, rho[1], rho[2], rho[1], 1, 0, rho[2], 0, 1), nrow=3)

  # Simulate random vector [G, H] using multivariate normal distribution
  sim.data <- mvtnorm::rmvnorm(N, mean = mean.vec, sigma = cov.mat)

  # Extract G and H from simulated data
  G <- sim.data[, 1]
  H <- sim.data[, -1]

  # Calculate indicator function for each simulation run
  indicators <- rep(0, N)
  for(i in 1:N){
    limits.list <- CI_limits(H[i,], b.knots.vec, b1.knots.vec,
                             s.knots.vec, r.knots.vec, t.vec, alpha, d, delta)
    if(limits.list$l <= G[i] & G[i] <= limits.list$u){
      indicators[i] <- 1
    }else{
      indicators[i] <- 0
    }

  }

  # Estimate coverage probability using Monte Carlo simulation
  cp.estimate <-  sum(indicators) / N

  # Standard error
  se <- sqrt(cp.estimate * (1 - cp.estimate) / N)

  # List to store the estimated coverage probability and it's standard error
  list(Estimate=cp.estimate, SE=se)
}
