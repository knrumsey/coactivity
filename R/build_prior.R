#' Build Prior Method for C_bass and Cfg_bass
#'
#' A method for building priors of the form needed in `C_bass`, `Cfg_bass` and similar functions.
#' For mixture distributions, all arguments (except for lower and upper) should be matrices (with nrow equal to the number of mixture components) rather than vectors.
#'
#'
#' @param dist A vector of length p. Valid entries include "uniform", "normal", "beta", "gamma".
#' @param lower A p-vector of lower truncation bounds. `-Inf` is a valid entry.
#' @param upper A p-vector of lower truncation bounds. `Inf` is a valid entry.
#' @param mu A p-vector of means (used for normal/truncated normal only)
#' @param sigma A p-vector of sds (used for normal/truncated normal only)
#' @param shape1 A p-vector of shape1 parameters for beta prior
#' @param shape2 A p-vector of shape2 parameters for beta prior
#' @param shape A p-vector of shape parameters for gamma prior
#' @param scale A p-vector of scale parameters for gamma prior
#' @param weights A vector of mixture weights of the same dimension as dist.
#' @return a list which can be passed into C_bass or Cfg_bass as a prior.
#' @details Builds a list for passing to the \code{coactivity} functions. List contains one component per input variable.
#' The \code{dist} argument must be passed in full, but all other values can be scalars (and will be reshaped accordingly).
#' Truncation bounds cannot vary by mixture component. See examples below.
#' @examples
#' # standard uniform priors for 5 inputs
#' build_prior(rep("uniform", 5), lower=0, upper=1)
#'
#' # truncated normals with different means for each input
#' mu_vec <- c(0.4, 0.5, 0.3, 0.7, 0.5)
#' build_prior(rep("normal", 5), lower=0, upper=1,
#'             mean=mu_vec, sd=0.1)
#'
#' # A mixture of normals (p=4)
#' mu_mat = matrix(c(0.4, 0.5, NA,
#'                 0.5, NA, NA,
#'                 0.25, 0.5, 0.75,
#'                 0.4, 0.5, 0.8),
#'                 ncol=3, byrow=TRUE)
#'weights_mat = matrix(c(2, 3, 0,
#'                 1, 0, 0,
#'                 1, 1, 1,
#'                 1, 1, 10),
#'                 ncol=3, byrow=TRUE)
#' build_prior(matrix("normal", nrow=4, ncol=3),
#'             lower=-Inf, upper=Inf,
#'             mean=mu_mat, sd=0.1,
#'             weights=weights_mat)
#'
#' @export
build_prior <- function(dist, lower=-Inf, upper=Inf,
                        mu=NULL, sigma=NULL,
                        shape1=NULL, shape2=NULL,
                        shape=NULL, scale=NULL,
                        weights=NULL) {
  # Determine dimensions
  if (is.null(dim(dist))) {
    p <- length(dist)
    K <- 1 # number of mixture components
    dist <- matrix(dist, nrow = p, ncol = K) # Expand dist to a matrix if not already
  } else {
    p <- nrow(dist)
    K <- ncol(dist)
  }

  # Helper function to recycle scalar values to appropriate size
  recycle_to_matrix <- function(arg, nrow, ncol) {
    if (is.null(arg)) return(NULL)
    if (is.vector(arg) && length(arg) == 1) {
      matrix(arg, nrow = nrow, ncol = ncol, byrow = TRUE)
    } else if (is.vector(arg) && length(arg) == ncol) {
      matrix(arg, nrow = nrow, ncol = ncol, byrow = FALSE)
    } else {
      arg
    }
  }

  # Expand scalar arguments as needed
  # Also expands
  lower <- recycle_to_matrix(lower, p, 1)
  upper <- recycle_to_matrix(upper, p, 1)
  mu <- recycle_to_matrix(mu, p, K)
  sigma <- recycle_to_matrix(sigma, p, K)
  shape1 <- recycle_to_matrix(shape1, p, K)
  shape2 <- recycle_to_matrix(shape2, p, K)
  shape <- recycle_to_matrix(shape, p, K)
  scale <- recycle_to_matrix(scale, p, K)
  weights <- recycle_to_matrix(weights, p, K)

  # Build the prior list
  prior <- list()
  for (i in 1:p) {
    distribution <- dist[i, , drop = FALSE]
    pr <- list(dist = distribution,
               trunc = c(lower[i], upper[i]))

    if (!is.null(mu)) {
      pr$mu <- mu[i, , drop = FALSE]
    }
    if (!is.null(sigma)) {
      pr$sigma <- sigma[i, , drop = FALSE]
    }
    if (!is.null(shape1)) {
      pr$shape1 <- shape1[i, , drop = FALSE]
    }
    if (!is.null(shape2)) {
      pr$shape2 <- shape2[i, , drop = FALSE]
    }
    if (!is.null(shape)) {
      pr$shape <- shape[i, , drop = FALSE]
    }
    if (!is.null(scale)) {
      pr$scale <- scale[i, , drop = FALSE]
    }
    if (!is.null(weights)) {
      weight_vec <- weights[i, , drop = FALSE]
      pr$weights <- weight_vec / sum(weight_vec)
    } else {
      pr$weights <- rep(1 / K, K) # Default equal weights if not specified
    }

    prior[[i]] <- pr
  }

  return(prior)
}
