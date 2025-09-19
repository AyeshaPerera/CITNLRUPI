#' Compute the coverage probability of the CIUUPI for a given vector of
#' \eqn{\gamma} and \eqn{\rho} using Theorem 4.5.1
#'
#' This function computes the coverage probability using Theorem 4.5.1 and Gauss
#' Legendre quadrature.
#' The coverage probability is
#' \deqn{
#' {\rm CP}(\boldsymbol{\gamma}, \boldsymbol{\rho}) = 1 - \alpha +
#' {\cal I}_1(\bm{\gamma}, \bm{\rho}) + {\cal I}_2(\bm{\gamma}, \bm{\rho})
#' }
#' We use Gauss Legendre quadrature to compute the integrals.
#'
#' @param gamma.vec Vector of \eqn{(\gamma_{p-1}, \gamma_{p})}
#' @param alpha The nominal coverage of the CI is \eqn{1 - \alpha}
#' @param d A positive number that specifies the interval \eqn{[-d,d]}
#' @param rho Vector of \eqn{(\rho_{p-1}, \rho_{p})}
#' @param t.vec The vector \eqn{(t_0, t_1, \dots, t_k)}
#' @param s.knots.vec The vector \eqn{(s(0), s(1), \dots, s(k))}
#' @param b.knots.vec The vector \eqn{(b(1), b(2), \dots, b(k))}
#' @param b1.knots.vec The vector \eqn{(b_1(0), b_1(1), \dots, b_1(k))}
#' @param r.knots.vec The vector \eqn{(r_1(0), r_1(1), \dots, r_1(k))}
#' @param delta A positive number that specifies the interval \eqn{[0,d]}
#' @param no.nodes.GL The number of Gauss Legendre nodes
#'
#' @returns Coverage probability computed using Theorem 4.5.1 and Gauss
#' Legendre quadrature.
#' @export
#'
#'
CP_T1_GL <- function(gamma.vec, alpha, d, rho, t.vec, s.knots.vec, b.knots.vec, b1.knots.vec, r.knots.vec, delta, no.nodes.GL){
  # This module computes the coverage probability
  #
  # Inputs
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies
  #    the interval [-d,d]
  # rho: vector of (rho{p-2}, rho{p})
  # t.vec: the vector (t0, t1, ..., tk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # b.knots.vec: the vector (b1, b2, ..., bk)
  # b1.knots.vec: the vector (b11, b12, ..., b1k)
  # r.knots.vec: the vector (r0, r1, ..., rk)
  #
  # Output
  # Coverage probability computed with integrate functions
  #
  # Written by A.Perera

  cp <- (1 - alpha) + ICP_T1_GL(rho, d, alpha, gamma.vec, t.vec,
           b.knots.vec, b1.knots.vec, s.knots.vec, r.knots.vec,
           delta, no.nodes.GL)

  cp

}
