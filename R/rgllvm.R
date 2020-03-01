#' rgllvm package
#'
#' This package provides functions to estimate parameters in a
#' generalized linear latent variable model where the latent variable
#' is a random intercept or random slope, and the distribution
#' of the latent variable is left unknown
#'
#' @references
#' Garcia, T.P. and Ma, Y. (2015). Optimal estimator for
#' logistic model with distribution-free random intercept.
#' Scandinavian Journal of Statistics, 43, 156-171.
#'
#' Ma, Y. & Genton, M. G. (2010). Explicit estimating equations for semiparametric generalized
#' linear latent variable models. Journal of the Royal Statistical Society, Series B 72, 475-495.
#'
#' Wei, Y., Ma, Y., Garcia, T.P. and Sinha, S. (2019). Consistent estimator
#' for logistic mixed effect models. The Canadian Journal of Statistics, 47, 140-156. doi:10.1002/cjs.11482.
#'
#' @useDynLib rgllvm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name rgllvm
NULL
