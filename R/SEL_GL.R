#' Compute the scaled expected length of the CIUUPI for a given vector of
#' \eqn{\gamma} using Gauss Legendre quadrature
#'
#' Evaluate the scaled expected length of the confidence interval that utilizes
#' uncertain prior information (CIUUPI) at a given \eqn{\gamma} using Gauss
#' Legendre quadrature. In the formula
#' \deqn{
#' {\rm SEL}(\bm{\gamma}) =
#' 1 + \frac{1}{ z_{1-\alpha/2}} \int_{-d}^d s(h_p)\,\phi(h_{p} - \gamma_{p})\,dh_p
#' + \frac{1}{ z_{1-\alpha/2}} \int_{-\delta}^\delta w_1(|h_p|)\,\phi(h_{p} - \gamma_{p})\,dh_p
#'   \int_{-d}^{d} r_1(h_{p-1})\,\phi(h_{p-1} - \gamma_{p-1})\,dh_{p-1}
#' - \big(\Phi(d - \gamma_p) - \Phi(-d - \gamma_p)\big).
#' }
#' all integrals are computed using Gauss Legendre quadrature.
#'
#' @param r.knots.vec The vector \eqn{(r_0, r_1, \dots, r_k)}
#' @param s.knots.vec The vector \eqn{(s_0, s_1, \dots, s_k)}
#' @param t.vec The vector \eqn{(t_0, t_1, \dots, t_k)}
#' @param alpha The nominal coverage of the CI is \eqn{1 - \alpha}
#' @param d A positive number that specifies the interval \eqn{[-d,d]}
#' @param gamma.vec Vector of \eqn{(\gamma_{p-1}, \gamma_{p})}
#' @param no.nodes.GL The number of Gauss Legendre nodes
#' @param delta A positive number that specifies the interval \eqn{[0,d]}
#'
#' @returns Scaled expected length computed using Gauss Legendre quadrature.
#' @export
#'
SEL_GL <- function(r.knots.vec, s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta){
  # This module computes the scaled expected length
  # using Gauss Legendre quadrature
  #
  # Inputs
  # r.knots.vec: the vector (r0, r1, ..., rk)
  # s.knots.vec: the vector (s0, s1, ..., sk)
  # t.vec: the vector (t0, t1, ..., tk)
  # alpha: The nominal coverage of the CI is 1 - alpha
  # d: a positive number that specifies
  #    the interval [-d,d]
  # gamma.vec: vector of (gamma{p-1}, gamma{p})
  # no.nodes.GL: the number of Gauss Legendre nodes
  #
  # Output
  # Scaled expected length
  #
  # Written by A.Perera August 2023

  1 + SEL_GL_t1(s.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL) + SEL_GL_t2(r.knots.vec, t.vec, alpha, d, gamma.vec, no.nodes.GL, delta) - (pnorm(d - gamma.vec[2]) - pnorm(-d - gamma.vec[2]))
}
