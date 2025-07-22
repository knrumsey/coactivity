#' Gradient evaluations of a BASS model
#'
#' This function behaves similarly to the predict.bass method, except that it returns
#' predictions for the gradient of f.
#'
#' @param object A fitted model, output from the bass function
#' @param newdata A matrix of input locations at which to evaluate the gradient. The columns should correspond to the same variables used in the \code{bass} function.
#' @param mcmc.use A vector indexing which MCMC iterations to be used.
#' @param verbose logical; should progress be displayed
#' @export
#' @return An array of gradient evaluations. First dimension indexes MCMC iteration, second dimension indexes input location, third dimension indexes the input variable.
#' @examples
#' #' # simulate data (Friedman function with first variable as functional)
#' n <- 500
#' p <- 2
#' X <- matrix(runif(n*p), ncol=p)
#' y <- apply(X, 1, duqling::dms_additive)
#' mod <- bass(X, y)
#'
#' Xnew <-  matrix(runif(100*p), ncol=p)
#' gradients <- gradient_bass(mod, Xnew)
#'
#' @export
gradient_bass <- function(object, newdata, mcmc.use=NULL, verbose=FALSE){
  if(is.null(mcmc.use)){
    mcmc.use <- seq_along(object$nbasis)
  }
  X <- BASS:::scale_range_mat(newdata, object$range.des)
  ranges <- apply(object$range.des, 2, diff)
  n <- nrow(X)
  p <- ncol(X)
  L <- length(mcmc.use)

  grads <- array(0, dim = c(L, n, p))
  for(ell in 1:L){
    if(verbose & (ell %% 100) == 0){
      cat("Starting MCMC iteration", ell, "\n")
    }
    mcmc_iter = mcmc.use[ell]
    model_number = object$model.lookup[mcmc_iter]

    # Number of basis functions
    M <- object$nbasis[mcmc_iter]
    beta <- object$beta[mcmc_iter, 2:(M+1)] # Skip intercept

    # Loop over each basis function
    for(m in 1:M){
      nint <- object$n.int.des[model_number, m]
      knotInds <- object$knotInd.des[model_number, m, 1:nint]
      signs <- object$signs.des[model_number, m, 1:nint]
      vars <- object$vars.des[model_number, m, 1:nint]
      knots <- sapply(seq_along(vars), function(i) object$xx.des[knotInds[i], vars[i]])

      # Account for Devin's scaling factors
      d <- 1/((signs + 1)/2 - signs*knots)
      signs <- signs*d

      # Loop only over variables in vars (u_{im} = 1)
      for(j in vars){
        ind <- which(j == vars)

        # Compute derivative only for variable j
        curr <- rep(NA, n)
        for(i in 1:n){
          # Compute h_jm'(x_j)

          if(abs(X[i, j] - knots[ind]) < 1e-4){
            browser()
          }

          h_prime <- beta[m] * signs[ind] * as.numeric(signs[ind] * (X[i, j] - knots[ind]) > 0)
          #h_prime <- beta[m] * signs[ind] * as.numeric(signs[ind] / abs(signs[ind]) * (X[i, j] - knots[ind]) > 0)

          # Compute product of h_km(x_k) for k != j
          vars_left <- vars[-ind]
          signs_left <- signs[-ind]
          knots_left <- knots[-ind]
          res_prod <- 1
          for(k in seq_along(vars_left)){
            jprime <- vars_left[k]
            #res_prod <- res_prod * as.numeric(signs_left[k] * (X[i, jprime] - knots_left[k]) > 0)
            hinge_val <- signs_left[k] * (X[i, jprime] - knots_left[k])
            res_prod  <- res_prod * pmax(hinge_val, 0)  # <-- keep full hinge
          }
          curr[i] <- h_prime * res_prod
        }
        grads[ell, , j] <- grads[ell, , j] + curr / ranges[j]
      }
    }
  }
  return(grads)
}
