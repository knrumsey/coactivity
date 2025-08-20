#' Convert functional BASS model to BASS model
#'
#' The argument to this function is the output of a \code{bass()} call when a single functional variable is specified using the \code{xx.func} argument.
#' Note that the resulting model may not be a valid bass object for some applications,
#' but the resulting model can be passed to \code{concordance::C_bass()} and related functions.
#'
#' @param bfm an object of class bass, where a functional variable has been specified.
#' @export
#' @examples
#' # simulate data (Friedman function with first variable as functional)
#' f<-function(x){
#' 10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#' }
#' sigma <- 1 # noise sd
#' n <- 500 # number of observations
#' nfunc <- 50 # size of functional variable grid
#' xfunc <- seq(0,1,length.out=nfunc) # functional grid
#' x <- matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
#' X <- cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
#' y <- matrix(f(X),nrow=n)+rnorm(n*nfunc,0,sigma)
#'
#' # fit BASS
#' mod_func <- bass(x,y,xx.func=xfunc)
#'
#' # convert to standard BASS model
#' mod <- bassfunc2bass(mod_func)
#'
#' # Estimate C matrix (augmentation approach)
#' C <- C_bass(mod)
#' @export
bassfunc2bass <- function(bfm){
  if(bfm$pfunc == 0){
    warning("No functional variable detected")
    return(bfm)
  }
  # First manipulate X and y into long form
  nx <- nrow(bfm$y)
  nt <- ncol(bfm$y)
  nmc <- length(bfm$nbasis)
  n <- length(bfm$y)
  pd <- bfm$pdes
  pf <- bfm$pfunc
  p <- pd + pf
  tt <- bfm$xx.func

  # Reformat the data
  XX <- cbind(
        kronecker(rep(1, nt), bfm$xx.des),
        kronecker(bfm$xx.func, rep(1, nx))
  )

  yy <- as.vector(bfm$y)
  if(!is.null(bfm$yhat)){
    yyhat.mean <- as.vector(bfm$yhat.mean)
    yyhat <- t(apply(bfm$yhat, 1, function(mat) as.vector(mat)))
  }else{
    yyhat <- NULL
    yyhat.mean <- NULL
  }


  # Get information about bass
  nm <- bfm$n.models
  nbasis <- bfm$nbasis
  maxBasis <- max(nbasis)
  maxInt <- bfm$maxInt.des + bfm$maxInt.func

  # Change the variable indicator for functional variable
  bfm$vars.func <- bfm$vars.func + pd

  # Combine information
  n.int <- bfm$n.int.des + bfm$n.int.func
  signs <- abind_custom(bfm$signs.des, bfm$signs.func, along=3)
  knots <- abind_custom(bfm$knotInd.des, bfm$knotInd.func, along=3)
  knotInds <- abind_custom(bfm$knotInd.des, bfm$knotInd.func, along=3)
  vars  <- abind_custom(bfm$vars.des, bfm$vars.func, along=3)
  range <- cbind(bfm$range.des, bfm$range.func)

  for(i in 1:nm){
    mcmc_iter <- which(bfm$model.lookup == i)
    nbasis_curr <- max(nbasis[mcmc_iter])
    for(j in seq_len(nbasis_curr)){
      # Fix knots so that they have data values, rather than indices
      knotInd_curr <- knotInds[i,j,]
      knot_curr <- knots[i,j,]
      vars_curr <- vars[i,j,]
      not_na_inds <- which(!is.na(knot_curr))
      if(maxInt %in% not_na_inds){
        # add time variable
        ind_old <- knot_curr[maxInt]
        ind_new <- (ind_old - 1) * nt + 1
        knot_curr[maxInt] <- tt[ind_old]
        knotInd_curr[maxInt] <- ind_new
        not_na_inds <- not_na_inds[-length(not_na_inds)]
      }
      if(length(not_na_inds) > 0){
        knot_inds <- cbind(knotInd_curr[not_na_inds], vars_curr[not_na_inds])
        knot_curr[not_na_inds] <- XX[knot_inds]
      }
      #if(any(knot_curr == 1, na.rm=TRUE)) browser()
      knots[i,j,] <- knot_curr

      # Shift vectors so that NAs occur at the end
      # This needs to be generalized if we want to allow for pfunc > 1
      nint_curr <- n.int[i,j]
      if(is.na(knot_curr[nint_curr])){
        knots[i,j,nint_curr] <- knots[i,j,maxInt]
        knots[i,j,maxInt]    <- NA

        knotInds[i,j,nint_curr] <- knotInds[i,j,maxInt]
        knotInds[i,j,maxInt]    <- NA

        signs[i,j,nint_curr] <- signs[i,j,maxInt]
        signs[i,j,maxInt]    <- NA

        vars[i,j,nint_curr] <- vars[i,j,maxInt]
        vars[i,j,maxInt]    <- NA
      }
    }
  }

  out <- bfm
  # Replace data fields
  out$xx.des <- XX
  out$y <- yy
  out$yhat <- yyhat
  out$yat.mean <- yyhat.mean
  out$p <- p

  # Add new fields
  out$n.int.des <- n.int
  out$signs.des <- signs
  out$vars.des  <- vars
  out$knots.des <- knots
  out$knotInd.des <- knotInds
  out$pdes      <- p
  out$maxInt.des<- maxInt
  out$range.des <- range
  out$cx <- c(out$cx, rep("numeric", pf))
  if(out$cat){
    out$type <- "_des_cat"
  }else{
    out$type <- "_des"
  }

  # Remove invalid fields
  out$des.basis <- NULL
  out$func.basis <- NULL
  out$curr.list <- NULL

  # Remove functional fields
  out$pfunc     <- 0
  out$n.int.func <- NULL
  out$signs.func <- NULL
  out$vars.func <- NULL
  out$knotInd.func <- NULL
  out$range.func <- NULL
  out$maxInt.func <- NULL
  out$func      <- FALSE
  out$wasfunc   <- TRUE # Use this to flag (in C_bass) that knots are not indices.
  class(out)    <- c("bass")
  return(out)
}
