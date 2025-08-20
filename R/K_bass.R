#' Estimate the Zahm Matrix with BASS
#'
#' Closed form estimator of the Zahm matrix using a BASS model
#'
#' @param mod a fitted functional BASS model, or a list of univariate bass models. The output of the bass() or bassPCA() functions.
#' @param prior a list, like one returned by the \code{build_prior()} function. See the documentation for details.
#' @param prior_func prior for the functional variable. Either a function over 0 to 1 (assumed to integrate to 1) or a vector with the same length as \code{func.use}.
#' @param mcmc.use a vector of indices telling which mcmc draws to use
#' @param use_native_scale logical (default `TRUE`). Determines the scale of the inputs for computing the \( C \)-matrix. When `TRUE`, the \( C \)-matrix is computed on the original (native) scale of the input variables. When `FALSE`, the \( C \)-matrix corresponds to the inputs normalized to the \([0, 1]\) range, as used internally by BASS. This also affects derived quantities, such as activity scores..
#' @param func.use a vector indicating which values of the functional variable to compute C for, if applicable
#' @param func.true An optional vector of values for the functional variable in bassPCA. Should have length equal to nrow(bassPCA$dat$basis).
#' @param pca_method takes value 1 or 2 to indicate which method should be used for estimating C. Ignored unless \code{class(mod) == "bassBasis"}.
#' @param verbose logical; should progress be displayed?
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details The K matrix is C(t) integrated over t with respect to a prior rho_t. In practice, this is equivalent to the Zahm matrix for multivariate response with a diagonal R.
#' @examples
#' # FRIEDMAN FUNCTION
#' # First input is treated as functional
#' # Use p=5, so there is one inert variable
#' f <- function(x, t) {
#'   10 * sin(pi * t * x[1]) + 20 * (x[2] - 0.5)^2 + 10 * x[3] + 5 * x[4]
#' }
#'
#' # ===========================================
#' #        GENERATE DATA
#' # ===========================================
#' XX <- lhs::randomLHS(500, 5)
#' y1 <- apply(XX, 1, f, t = 0.5)
#' xfunc <- seq(0, 1, length.out = 20)
#' yfunc <- t(apply(XX, 1, f, t = xfunc))
#'
#' Xtest <- lhs::randomLHS(100, 5)
#' ytest <- apply(Xtest, 1, f, t = 0.5)
#'
#' # ===========================================
#' #        CASE 1: Univariate BASS
#' # ===========================================
#' mod1 <- bass(XX, y1)
#' C1 <- C_bass(mod1)
#'
#' # ===========================================
#' #      CASE 2: Augmented BASS (fixed t)
#' # ===========================================
#' mod2_full <- bass(XX, yfunc, xx.func = xfunc)
#' mod2 <- bassfunc2bass_fixed_t(mod2_full, 0.5)
#' C2 <- C_bass(mod2)
#' C2b <- C_bass(mod2_full, func.use=0.5)
#'
#' # ===========================================
#' #      CASE 3: PCA BASS (fixed t)
#' # ===========================================
#' mod3_full <- bassPCA(XX, yfunc)
#' mod3 <- bassPCA2bass_fixed_t(mod3_full, 0.5)
#' C3 <- C_bass(mod3)
#' C3b <- C_bass(mod3_full, func.use=0.5)
#'
#' @export
K_bass <- function(mod, prior=NULL, prior_func=function(tt) dunif(tt),
                   mcmc.use=NULL, use_native_scale=FALSE, func.use=NULL,
                   func.true=NULL, pca_method=2,
                   verbose=FALSE){
  # Determine what to do with func.use
  if(is.null(func.use)){
    if(class(mod) == "list"){
      nfunc <- length(mod)
      func.use <- seq(0, 1, length.out=nfunc)
    }
    if(class(mod) == "bassBasis"){
      nfunc <- nrow(mod$dat$basis)
      func.use <- seq(0, 1, length.out=nfunc)
    }
    if("bass" %in% class(mod)){
      if(mod$func){
        nfunc <- nrow(mod$xx.func)
        func.use <- as.numeric(mod$xx.func)
      }
    }
  }else{
    nfunc <- length(func.use)
    if(nfunc == 1){
      warning("K_bass doesn't make sense with just one time point.")
    }
  }

  # Handle prior
  if(is.function(prior_func)){
    rho_t <- prior_func(func.use)
  }else{
    rho_t <- prior_func
  }
  if(length(rho_t) != nfunc){
    stop("Prior for functional variables is the wrong size.")
  }
  rho_t <- rho_t/sum(rho_t)

  # Calculate C Matrices
  res <- C_bass(mod, prior, mcmc.use, use_native_scale, func.use, func.true, pca_method, verbose)

  if(is.list(res)){
    # Combine with prior weights
    if(length(res) != length(rho_t)){
      stop("prior length doesn't match length of output")
    }
  }else{
    if(length(rho_t) != 1){
      stop("prior length doesn't match length of output")
    }
  }

  # Get nmc and dimension of matrix
  K <- list()
  if(nfunc > 1){
    if(is.list(res[[1]])){
      nmc <- length(res[[1]])
    }else{
      nmc <- 1
    }
    if(is.list(res[[1]])){
      p <- nrow(res[[1]][[1]])
    }else{
     p <- nrow(res[[1]])
    }
  }else{
    # nfunc == 1, bizarre case.
    # Lets just fill it out and return here
    if(is.list(res)){
      nmc <- length(res)
    }else{
      nmc <- 1
    }
    if(is.list(res)){
      p <- nrow(res[[1]])
      K <- list()
      for(j in 1:nmc){
        K[[j]] <- rho_t * res[[j]]
      }
      if(nmc == 1){
        return(K[[1]])
      }else{
        return(K)
      }
    }else{
      p <- nrow(res)
      K <- rho_t * res
      return(K)
    }
  }

  # Initialize K
  K <- list()
  for(j in seq_len(nmc)){
    K[[j]] <- matrix(0, nrow=p, ncol=p)
  }

  # Fill out K
  # We can assume nfunc > 1
  for(i in seq_len(nfunc)){
    for(j in seq_len(nmc)){
      curr <- res[[i]]
      if(nmc > 1){
        curr <- curr[[j]]
      }
      K[[j]] <- K[[j]] + rho_t[i] * curr
    }
  }

  if(nmc == 1){
    return(K[[1]])
  }else{
    return(K)
  }
}



