#' Compute the coverage probability of the CIUUPI for a given vector of
#' \eqn{\gamma} and \eqn{\rho} using Theorem 4.5.2
#'
#' This function computes the coverage probability using Theorem 4.5.2 and
#' multiple \code{integrate()} functions.
#' The coverage probability is
#' \deqn{
#' {\rm CP}(\boldsymbol{\gamma}, \boldsymbol{\rho}) = 1 - \alpha - t_1(\gamma_p)
#' + t_2(\bm{\gamma}; b, s, b_1, r_1) + t_3(\bm{\gamma}; b, s)
#' + t_4(\gamma_p; b, s)
#' }
#' We use the \code{integrate()} function to compute the integrals.
#'
#' @param rho Vector of \eqn{(\rho_{p-1}, \rho_{p})}
#' @param d A positive number that specifies the interval \eqn{[-d,d]}
#' @param alpha The nominal coverage of the CI is \eqn{1 - \alpha}
#' @param gamma.vec Vector of \eqn{(\gamma_{p-1}, \gamma_{p})}
#' @param t.vec The vector \eqn{(t_0, t_1, \dots, t_k)}
#' @param b.knots.vec The vector \eqn{(b(1), b(2), \dots, b(k))}
#' @param b1.knots.vec The vector \eqn{(b_1(0), b_1(1), \dots, b_1(k))}
#' @param s.knots.vec The vector \eqn{(s(0), s(1), \dots, s(k))}
#' @param r.knots.vec The vector \eqn{(r_1(0), r_1(1), \dots, r_1(k))}
#' @param delta A positive number that specifies the interval \eqn{[0,d]}
#'
#' @returns Coverage probability computed using Theorem 4.5.2 and the
#' \code{integrate()} function.
#' @export
#'
#'
CP_T2_int <- function(rho, d, alpha, gamma.vec, t.vec, b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec, delta){
  # Computes coverage probability using integrate function
  #
  # Input
  # rho: vector of (rho{p-2}, rho{p})
  # d: a positive number that specifies the interval [-d,d]
  # alpha: The nominal coverage of the CI is 1 - alpha
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # t.vec: the vector (t0, t1, ..., tk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # delta: tuning constant for w1(x) were 0 < delta <= d
  # no.nodes.GL
  #
  # Output
  # Coverage probability
  #
  # Written by A. Perera June 2023

  z <- qnorm(1 - alpha / 2)
  mu <- c(0, gamma.vec[2])
  Sigma <- matrix(c(1, rho[2], rho[2], 1), nrow=2)

  term1 <- pmvnorm(lower=c(-z, -d), upper=c(z, d),
           algorithm = "Miwa", mean=mu, sigma=Sigma)[1]

  term2 <- k_phi_int_i(rho, d, alpha, gamma.vec, t.vec,
           b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec,
           delta)

  term3 <- CP_T2_int_t3(rho, gamma.vec, alpha, d, t.vec,
           b.knots.vec, s.knots.vec, delta)

  term4 <- CP_T2_int_t4(rho, d, alpha, gamma.vec, t.vec,
           b.knots.vec, s.knots.vec, delta)

  1 - alpha - term1 + term2 + term3 + term4

}
