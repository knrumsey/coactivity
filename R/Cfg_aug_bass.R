#' Estimate the Augmented Co-Constantine Matrix
#'
#' Closed form estimator of the Cfg_aug matrix using a functional BASS model
#'
#' @param mod1 first fitted (functional) BASS model. The output of the bass() or bassPCA() functions.
#' @param mod2 second fitted (functional) BASS model. The output of the bass() or bassPCA() functions.
#' @param prior a list, like one returned by the \code{build_prior()} function. See the documentation for details.
#' @param mcmc.use a vector of indices telling which mcmc draws to use
#' @param use_native_scale logical (default `TRUE`). Determines the scale of the inputs for computing the \( C \)-matrix. When `TRUE`, the \( C \)-matrix is computed on the original (native) scale of the input variables. When `FALSE`, the \( C \)-matrix corresponds to the inputs normalized to the \([0, 1]\) range, as used internally by BASS. This also affects derived quantities, such as activity scores..
#' @param func.true What are the true values of the functional variable, for a bassPCA model?
#' @param verbose logical; should progress be displayed?
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details The posterior distribution of the augmented C matrix. Note that mod must be a functional bass model; either a \code{bass()} call with \code{xx.func} specified, or an object of class \code{bassBasis} such as the output of a \code{bassPCA()} call.
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
#' @export
Cfg_aug_bass <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.true=NULL, verbose=FALSE){
  if("bass" %in% class(mod1)){
    if(!("bass" %in% class(mod2))){
      stop("mod1 and mod2 must have the same class")
    }
    if(mod1$pfunc == 0 | mod2$pfunc == 0){
      stop("models must be functional BASS models.")
    }
    mod_new1 <- bassfunc2bass(mod1)
    mod_new2 <- bassfunc2bass(mod2)
    res <- Cfg_bass(mod_new1, mod_new2, prior, mcmc.use, use_native_scale, verbose=verbose)
  }else if(class(mod1) == "bassBasis"){
    if(class(mod2) != "bassBasis"){
      stop("mod1 and mod2 must have the same class")
    }
    nrow1 <- nrow(mod1$dat$basis)
    nrow2 <- nrow(mod2$dat$basis)
    if(nrow1 != nrow2){
      if(is.null(func.true)){
        warning("num rows of mod1 and mod2 don't match. Consider specifying func.true to avoid ambiguity.")
      }

    }
    nfunc <- nrow1
    if(is.null(func.true)){
      func.true <- seq(0, 1, length.out=nfunc)
    }
    if(is.null(prior)){
      rho_t <- rep(1, nfunc)/nfunc
    }else{
      p <- length(prior) - 1
      pf <- prior[[p+1]]
      if(pf$dist == "uniform"){
        rho_t <- dunif(func.true, pf$trunc[0], pf$trunc[1])
        scale_t <- 1
      }else if(pf$dist == "beta"){
        rho_t <- dbeta(func.true, pf$shape1, pf$shape2)
        scale_t <- 1
      }else{
        warning("prior for functional variable should be on zero to 1.")
        scale_t <- diff(range(func.true))
        if(pf$dist == "normal"){
          rho_t <- dnorm(func.true, pf$mu, pf$sigma)
        }else if(pf$dist == "gamma"){
          rho_t <- dgamma(func.true, pf$shape, scale=pf$scale)
        }else{
          warning("prior for t not recognized. Using uniform")
          rho_t <- dunif(func.true, 0, 1)
        }
      }
    }
    res <- Cfg_aug_bassPCA(mod1, mod2, prior, rho_t, mcmc.use, use_native_scale, func.true, verbose)
  }else{
    stop("mod must be a functional BASS model")
  }
  return(res)
}

Cfg_aug_bassPCA <- function(mod1, mod2, prior, rho_t, mcmc.use, use_native_scale, func.true, verbose){
  mod <- mod1
  if(is.null(mcmc.use)){
    mcmc.use <- length(mod$mod.list[[1]]$nbasis)
  }
  nmc <- length(mcmc.use)

  # Model 1
  U1 <- mod$dat$basis
  mod_list1 <- mod$mod.list
  K1 <- length(mod_list1)
  p <- mod_list1[[1]]$pdes

  # Model 2
  U2 <- mod2$dat$basis
  mod_list2 <- mod2$mod.list
  K2 <- length(mod_list2)
  p2 <- mod_list2[[1]]$pdes
  if(p != p2) warning("number of variables doesn't match in mod1 and mod2")

  if(nmc == 1){
    Caug <- matrix(0, nrow=p+1, ncol=p+1)
  }else{
    Caug <- list()
    for(i in seq_len(nmc)){
      Caug[[i]] <- matrix(0, nrow=p+1, ncol=p+1)
    }
  }
  for(k1 in 1:K1){
    for(k2 in 1:K2){
      obj <- Cfg_aug_bass_univariate(mod_list1[[k1]], mod_list2[[k2]], prior, mcmc.use, verbose)

      u1 <- U1[,k1]
      u2 <- U2[,k2]
      u1_deriv <- compute_derivative(u1, func.true)
      u2_deriv <- compute_derivative(u2, func.true)
      ups00 <- trapezoidal(func.true, u1 * u2 * rho_t)
      ups01 <- trapezoidal(func.true, u1_deriv * u2 * rho_t)
      ups10 <- trapezoidal(func.true, u1 * u2_deriv * rho_t)
      ups11 <- trapezoidal(func.true, u1_deriv * u2_deriv * rho_t)

      Cfg <- list_multiply(ups00, obj$Cfg)
      Lfg <- list_multiply(ups01, obj$Lfg)
      Lgf <- list_multiply(ups10, obj$Lgf)
      lfg <- list_multiply(ups11, obj$lfg)

      if(nmc == 1){
        Caug_curr <- cbind(rbind(Cfg, Lfg), c(Lgf, lfg))
        Caug <- Caug + Caug_curr
      }else{
        for(i in 1:nmc){
          Caug_curr <- cbind(rbind(Cfg[[i]], Lfg[[i]]), c(Lgf[[i]], lfg[[i]]))
          Caug[[i]] <- Caug[[i]] + Caug_curr
        }
      }
    }
  }
  return(Caug)
}
