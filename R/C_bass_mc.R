#' Estimate the Constantine Matrix with BASS using Monte Carlo
#'
#' Monte Carlo Estimator of the C matrix using a BASS model
#'
#' @param mod a fitted BASS model. The output of the bass() or bassPCA() functions.
#' @param prior A function providing samples of x. Or, alternatively, a list, like one returned by the \code{build_prior()} function. See the documentation for details.
#' @param n_mc number of Monte Carlo samples
#' @param mcmc.use a vector of indices telling which mcmc draws to use
#' @param use_native_scale logical (default `TRUE`). Determines the scale of the inputs for computing the \( C \)-matrix. When `TRUE`, the \( C \)-matrix is computed on the original (native) scale of the input variables. When `FALSE`, the \( C \)-matrix corresponds to the inputs normalized to the \([0, 1]\) range, as used internally by BASS. This also affects derived quantities, such as activity scores..
#' @param func.use a vector indicating which values of the functional variable to compute C for, if applicable
#' @param func.true An optional vector of values for the functional variable in bassPCA. Should have length equal to nrow(bassPCA$dat$basis).
#' @param pca_method takes value 1 or 2 to indicate which method should be used for estimating C. See details.
#' @param verbose logical; should progress be displayed?
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details The C matrices are computed using inputs which are scaled to (0, 1). The \code{use_native_scale} flag indicates whether the C matrix should be transformed to the native space before returning.
#'
#' The \code{func.use} argument is used only if \code{mod} is fit to functional response data (using \code{bass()} with \code{xx.func} specified or by using \code{bassPCA}). In the latter case, the functional input is assumed to be between 0 and 1.
#'
#' The \code{pca_method} argument is used only when \code{mod} has class bassBasis. When set to 1 (the default), the decomposition theorem is used to estimate C. Otherwise, the model list is converted to a single model using \code{basslc2bass} and C is found directly. Method 1 is typically faster (especially for multiple func.use), but method 2 gives more flexibility in the choice of prior.
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
#' C1 <- C_bass_mc(mod1)
#'
#' @export
C_bass_mc <- function(mod, n_mc = 1e4, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.use=NULL, func.true=NULL, pca_method=2, verbose=FALSE){
  if(class(mod) == "list"){
    nmods <- length(mod)
    res <- list()
    for(i in seq_len(nmods)){
      res[[i]] <- C_bass_mc(mod[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose, n_mc = 1e4)
    }
    return(res)
  }
  if("bass" %in% class(mod)){
    if(mod$func){
      if(is.null(func.use)){
        func.use <- mod$xx.func
      }
      # Convert to univariate bass models
      mod_t <- bassfunc2bass_fixed_t(mod, func.use, verbose)
      if(length(func.use) > 1){
        res <- list()
        for(i in seq_along(func.use)){
          res[[i]] <- C_bass_univariate_mc(mod_t[[i]], prior, mcmc.use, use_native_scale, verbose, n_mc = 1e4)
        }
      }else{
        res <- C_bass_univariate_mc(mod_t, prior, mcmc.use, use_native_scale, verbose, n_mc = 1e4)
      }
      return(res)
    }else{
      res <- C_bass_univariate_mc(mod, prior, mcmc.use, use_native_scale, verbose, n_mc = 1e4)
      return(res)
    }
  }else if("bassBasis" %in% class(mod)){
    if(pca_method == 1){
      C_bass_pca_v1_mc(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose, n_mc = 1e4)
    }else{
      C_bass_pca_v2_mc(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose, n_mc = 1e4)
    }

  }else{
    return("mod should be a fitted bass model")
  }
}

C_bass_pca_v1_mc <- function(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose, n_mc = 1e4){
  warning("Not currently supported. Use pca_method=2")
  return(TRUE)
}

C_bass_pca_v2_mc <- function(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose, n_mc = 1e4){
  mod_list <- mod$mod.list
  K <- mod$dat$n.pc
  U <- mod$dat$basis
  nt <- nrow(U)
  # Convert models
  mod_u_list <- bassPCA2bass_fixed_t(mod, func.use, func.true)

  res <- list()
  cnt <- 1
  for(t_val in func.use){
    # Convert linear combination to bass
    mod_curr <- mod_u_list[[cnt]]
    res[[cnt]] <- C_bass_univariate_mc(mod_curr, prior, mcmc.use, use_native_scale, verbose, n_mc = 1e4)
    cnt <- cnt + 1
  }

  if(length(func.use) == 1){
    return(res[[1]])
  }else{
    return(res)
  }
}

C_bass_univariate_mc <- function(mod, prior = NULL, mcmc.use=NULL, use_native_scale=FALSE, verbose=TRUE, n_mc = 1e4){
  if(!is.null(mod$pfunc)){
    if(mod$pfunc > 0 && !isTRUE(mod$wasfunc)){
      warning("A functional variable was detected. Model will be converted. Use function bassfunc2bass to surpress this warning.")
      mod <- bassfunc2bass(mod)
    }
  }
  if(is.null(mcmc.use)){
    mcmc.use <- length(mod$nbasis)
  }
  gbass_flag <- "gbass" %in% class(mod)
  bassfunc_flag <- isTRUE(mod$wasfunc)
  if(gbass_flag || bassfunc_flag){
    mod$knotInd.des <- mod$knots.des
  }else{
    Xt <- mod$xx.des
  }

  # Parse the prior information
  if(is.null(prior)){
    prior <- list()
    for(ii in 1:(mod$pdes)){
      prior[[ii]] <- list(dist="uniform")
    }
  }
  if(is.list(prior)){
    #HANDLE PRIORS
    for(i in 1:(mod$pdes)){
      # 1. Scale prior truncation bounds to BASS scale
      if(is.null(prior[[i]]$trunc)){
        prior[[i]]$trunc <- mod$range.des[,i]
      }
      prior[[i]]$trunc <- scale_range(prior[[i]]$trunc, mod$range.des[,i])

      # 2. Handle each distribution type separately
      distribution = prior[[i]]$dist
      if(distribution == "normal"){
        if(is.null(prior[[i]]$weights)){
          num_mix <- length(prior[[i]]$mean)
          prior[[i]]$weights <- rep(1/num_mix, num_mix)
        }
        prior[[i]]$mean <- scale_range(prior[[i]]$mean, mod$range.des[,i])
        prior[[i]]$sd <- prior[[i]]$sd/(mod$range.des[2,i] - mod$range.des[1,i])
        prior[[i]]$z <- pnorm((prior[[i]]$trunc[2]-prior[[i]]$mean)/prior[[i]]$sd) - pnorm((prior[[i]]$trunc[1]-prior[[i]]$mean)/prior[[i]]$sd)
        cc <- sum(prior[[i]]$weights*prior[[i]]$z)
        prior[[i]]$weights <- prior[[i]]$weights/cc
        # DF: prior[[i]]$z
        # change weights with truncation
        # divide by cc instead to keep the same prior shape
        # does the truncation change the distribution shape in the non-truncated regions??

        #Check for extrapolation
        qq <- qnorm(c(0.001, 0.999), prior[[i]]$mean, prior[[i]]$sd)
        if((qq[1] < 0 - 0.1) | (qq[2] > 1 + 0.1)){
          warning('You are asking the emulator to extrapolate. This is not reccomended.')
        }
      }
      if(distribution == "uniform"){
        prior[[i]]$weights = 1
        #Check for extrapolation
        if(((prior[[i]]$trunc[1] < 0 - 0.1) | (prior[[i]]$trunc[2] > 1 + 0.1))){
          warning('You are asking the emulator to extrapolate. This is not reccomended.')
        }
      }
      if(distribution == "beta"){
        prior[[i]]$weights = 1
        if(abs(min(mod$xx.des[,i]) - 0) > 0.01 | abs(max(mod$xx.des[,i]) - 1) > 0.01){
          warning("For dist=\"beta\", data should be scaled to (0, 1) range BEFORE fitting bass model")
        }
      }
      if(distribution == "gamma"){
        prior[[i]]$weights = 1
        if(abs(min(mod$xx.des[,i]) - 0) > 0.01){
          warning("For dist=\"gamma\", data should be shifted to have a minimum of 0 BEFORE fitting bass model")
        }
        prior[[i]]$rate = prior[[i]]$rate * prior[[i]]$trunc[2]
      }
    }

    #browser()
    # Make prior into a function
    rprior <- get_prior_generator(prior)
  }else{
    rprior <- prior
  }

  # Start monte carlo loop
  X_mc <- t(sapply(1:n_mc, function(xx) rprior()))
  grads <- gradient_bass(mod, X_mc, mcmc.use=mcmc.use)


  n_mcmc <- dim(grads)[1]
  p <- dim(grads)[3]
  C_out <- list()
  A_tform <- diag(1/apply(mod$range.des, 2, diff))
  for(i in seq_along(mcmc.use)){
    C_tmp <- matrix(0, nrow=p, ncol=p)
    for(j in seq_len(n_mc)){
      C_tmp <- C_tmp + tcrossprod(grads[i,j,]) / n_mc
    }
    if(use_native_scale){
      # Transform back to native space
      C_tmp <- t(A_tform)%*%C_tmp%*%A_tform
    }
    C_out[[i]] <- C_tmp
  }
  if(n_mcmc == 1){
    return(C_out[[1]])
  }
  return(C_out)
}

#' Estimate the Constantine Matrix with BASS using Monte Carlo
#'
#' Monte Carlo estimator of the C matrix using a BASS model
#'
#' @param mod1 a fitted BASS model for the first function.
#' @param mod2 a fitted BASS model for the second function.
#' @param prior a list, like one returned by the \code{build_prior()} function. See the documentation for details.
#' @param mcmc.use a two-column matrix of MCMC indices corresponding to mod1 and mod2 respectively.
#' @param use_native_scale logical (default `TRUE`). Determines the scale of the inputs for computing the \( C \)-matrix. When `TRUE`, the \( C \)-matrix is computed on the original (native) scale of the input variables. When `FALSE`, the \( C \)-matrix corresponds to the inputs normalized to the \([0, 1]\) range, as used internally by BASS. This also affects derived quantities, such as activity scores..
#' @param func.use a vector indicating which values of the functional variable to compute C for, if applicable
#' @param pca_method takes value 1 or 2 to indicate which method should be used for estimating C. See details.
#' @param verbose logical; should progress be displayed?
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details The C matrices are computed using inputs which are scaled to (0, 1). The \code{use_native_scale} flag indicates whether the C matrix should be transformed to the native space before returning.
#'
#' The \code{func.use} argument is used only if \code{mod} is fit to functional response data (using \code{bass()} with \code{xx.func} specified or by using \code{bassPCA}). In the latter case, the functional input is assumed to be between 0 and 1.
#'
#' The \code{pca_method} argument is used only when \code{mod} has class bassBasis. When set to 1 (the default), the decomposition theorem is used to estimate C. Otherwise, the model list is converted to a single model using \code{lcbass2bass} and C is found directly. Method 1 is typically faster (especially for multiple func.use), but method 2 gives more flexibility in the choice of prior.
#' @examples
#' # FRIEDMAN FUNCTION AND SECONDARY FUNCTION
#' # First input is treated as functional
#' # Use p=5, so there is one inert variable
#' f <- function(x, t) {
#'   10 * sin(pi * t * x[1]) + 20 * (x[2] - 0.5)^2 + 10 * x[3] + 5 * x[4]
#' }
#' g <- function(x, t) {
#'   5 * sin(2 * pi * t * x[1]) + 15 * (x[3] - 0.3)^2 + 8 * x[4] + 3 * x[5]
#' }
#'
#' # ===========================================
#' #        GENERATE DATA
#' # ===========================================
#' XX <- lhs::randomLHS(500, 5)
#' y1_f <- apply(XX, 1, f, t = 0.5)
#' y1_g <- apply(XX, 1, g, t = 0.5)
#' xfunc <- seq(0, 1, length.out = 20)
#' yfunc_f <- t(apply(XX, 1, f, t = xfunc))
#' yfunc_g <- t(apply(XX, 1, g, t = xfunc))
#'
#' Xtest <- lhs::randomLHS(100, 5)
#' ytest_f <- apply(Xtest, 1, f, t = 0.5)
#' ytest_g <- apply(Xtest, 1, g, t = 0.5)
#'
#' # ===========================================
#' #        CASE 1: Univariate BASS
#' # ===========================================
#' mod1_f <- bass(XX, y1_f)
#' mod1_g <- bass(XX, y1_g)
#' C1_fg <- Cfg_bass_mc(mod1_f, mod1_g)
#' @export
Cfg_bass_mc <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.use=NULL, func.true=NULL, pca_method=2, verbose=FALSE){
  # HANDLE CASE WHERE mod1 or mod2 IS A LIST
  if(class(mod1) == "list" && class(mod2) == "list"){
    nmods <- length(mod1)
    res <- list()
    for(i in seq_len(nmods)){
      res[[i]] <- Cfg_bass_mc(mod1[[i]], mod2[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
    }
    return(res)
  }else{
    if(class(mod1) == "list"){
      nmods <- length(mod1)
      res <- list()
      for(i in seq_len(nmods)){
        res[[i]] <- Cfg_bass_mc(mod1[[i]], mod2, prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
      }
      return(res)
    }
    if(class(mod1) == "list"){
      nmods <- length(mod2)
      res <- list()
      for(i in seq_len(nmods)){
        res[[i]] <- Cfg_bass_mc(mod1, mod2[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
      }
      return(res)
    }
  }

  if (!("bass" %in% class(mod1) || "bassBasis" %in% class(mod1))) {
    stop("mod1 must be either a fitted BASS model or a BASS PCA model.")
  }
  if (!("bass" %in% class(mod2) || "bassBasis" %in% class(mod2))) {
    stop("mod2 must be either a fitted BASS model or a BASS PCA model.")
  }

  # Case 1: Both models are scalar BASS (pfunc = 0)
  if ("bass" %in% class(mod1) && "bass" %in% class(mod2) && !mod1$func && !mod2$func) {
    return(Cfg_bass_univariate_mc(mod1, mod2, prior, mcmc.use, use_native_scale, verbose))
  }

  # Case 2: Both models are functional BASS (pfunc > 0)
  if ("bass" %in% class(mod1) && "bass" %in% class(mod2) && mod1$func && mod2$func) {
    if (is.null(func.use)) func.use <- mod1$xx.func
    nfunc <- length(func.use)
    res <- list()
    cnt <- 1
    for (t_val in func.use) {
      mod1_t <- bassfunc2bass_fixed_t(mod1, t_val, verbose)
      mod2_t <- bassfunc2bass_fixed_t(mod2, t_val, verbose)
      res[[cnt]] <- Cfg_bass_univariate_mc(mod1_t, mod2_t, prior, mcmc.use, use_native_scale, verbose)
      cnt <- cnt + 1
    }
    if(nfunc == 1) res <- res[[1]]
    return(res)
  }

  # Case 3: One functional BASS (pfunc > 0) and one scalar BASS (pfunc = 0)
  if ("bass" %in% class(mod1) && "bass" %in% class(mod2) && xor(mod1$func, mod2$func)) {
    if(mod1$func){
      func_mod <- mod1
      scal_mod <- mod2
      swap_ord <- FALSE
    }else{
      func_mod <- mod2
      scal_mod <- mod1
      swap_ord <- TRUE
    }
    if (is.null(func.use)) func.use <- func_mod$xx.func
    nfunc <- length(func.use)
    res <- list()
    cnt <- 1
    for (t_val in func.use) {
      if(verbose & nfunc > 5 & cnt %in% round(seq(0, nfunc, length.out=5))[-1]){
        cat("Processing func value ", t_val, ". ", myTimestamp(), "\n", sep="")
      }
      func_mod_t <- bassfunc2bass_fixed_t(func_mod, t_val, verbose)
      if(swap_ord){
        res[[cnt]] <- Cfg_bass_univariate_mc(scal_mod, func_mod_t, prior, mcmc.use, use_native_scale, verbose)
      }else{
        res[[cnt]] <- Cfg_bass_univariate_mc(func_mod_t, scal_mod, prior, mcmc.use, use_native_scale, verbose)
      }
      cnt <- cnt + 1
    }
    if(nfunc == 1) res <- res[[1]]
    return(res)
  }

  # Case 4: At least one of the models is BASS PCA
  if ("bassBasis" %in% class(mod1) || "bassBasis" %in% class(mod2)) {
    if (is.null(func.use)) func.use <- mod1$xx.func
    if(pca_method == 1){
      return(Cfg_bass_pca_v1_mc(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose))
    }else{
      return(Cfg_bass_pca_v2_mc(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose))
    }
  }
}

Cfg_bass_pca_v1_mc <- function(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
  warning("Not currently supported. Please use pca_version = 2.")
  return(TRUE)
}

Cfg_bass_pca_v2_mc <- function(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
  nfunc <- length(func.use)
  if("bassBasis" %in% class(mod1)){
    type1 <- "pca"
  }else if("bass" %in% class(mod1)){
    if(mod1$func){
      type1 <- "func"
    }else{
      type1 <- "bass"
    }
  }else{
    stop("mod1 is not a valid argument")
  }

  if("bassBasis" %in% class(mod2)){
    type2 <- "pca"
  }else if("bass" %in% class(mod2)){
    if(mod2$func){
      type2 <- "func"
    }else{
      type2 <- "bass"
    }
  }else{
    stop("mod1 is not a valid argument")
  }

  # Make list of models for each t_val
  if(type1 == "pca"){
    model_list1 = bassPCA2bass_fixed_t(mod1, func.use=func.use, func.true=func.true)
  }
  if(type1 == "func"){
    model_list1 = bassfunc2bass_fixed_t(mod1, func.use=func.use)
  }
  if(type2 == "pca"){
    model_list2 = bassPCA2bass_fixed_t(mod2, func.use=func.use, func.true=func.true)
  }
  if(type2 == "func"){
    model_list2 = bassfunc2bass_fixed_t(mod2, func.use=func.use)
  }

  res <- list()
  cnt <- 1
  for(t_val in func.use){
    if(type1 == "bass"){
      curr1 <- mod1
    }else{
      if(nfunc == 1){
        curr1 <- model_list1
      }else{
        curr1 <- model_list1[[cnt]]
      }
    }
    if(type2 == "bass"){
      curr2 <- mod2
    }else{
      if(nfunc == 1){
        curr2 <- model_list2
      }else{
        curr2 <- model_list2[[cnt]]
      }
    }

    res[[cnt]] <- Cfg_bass_univariate_mc(curr1, curr2, prior, mcmc.use, use_native_scale, verbose)
    cnt <- cnt + 1
  }

  if(nfunc == 1){
    return(res[[1]])
  }else{
    return(res)
  }
}

Cfg_bass_univariate_mc <- function(mod1, mod2, prior = NULL, mcmc.use=NULL, use_native_scale=FALSE, verbose=TRUE, n_mc = 1e4){
  mod <- mod1
  if(!is.null(mod$pfunc)){
    if(mod$pfunc > 0 && !mod$wasfunc){
      warning("A functional variable was detected. Models will be converted. Use function bassfunc2bass to surpress this warning.")
      mod <- bassfunc2bass(mod)
      mod2 <- bassfunc2bass(mod2)
    }
  }

  # Check that the covariate matrix is the same
  if(max(abs(mod1$xx.des - mod2$xx.des)) > 1e-9){
    warning("Covariates are different (or are ordered differently) for mod1 and mod2. Results may be unreliable.")
  }
  if(is.null(mcmc.use)){
    mcmc.use <- min(length(mod$s2), length(mod2$s2))
  }
  mcmc.use <- as.matrix(mcmc.use)
  if(ncol(mcmc.use) == 1){
    mcmc.use <- cbind(mcmc.use, mcmc.use)
  }
  if(ncol(mcmc.use) > 2){
    warning("ncol(mcmc.use) should not exceed 2")
  }
  if(mod$pdes != mod2$pdes){
    stop("Detected different number of variables in mod1 and mod2")
  }
  gbass_flag1 <- "gbass" %in% class(mod)
  bassfunc_flag1 <- isTRUE(mod$wasfunc)
  if(gbass_flag1 || bassfunc_flag1){
    mod$knotInd.des <- mod$knots.des
  }

  gbass_flag2 <- "gbass" %in% class(mod2)
  bassfunc_flag2 <- isTRUE(mod2$wasfunc)
  if(gbass_flag2 || bassfunc_flag2){
    mod2$knotInd.des <- mod2$knots.des
  }
  # Parse the prior information
  if(is.null(prior)){
    prior <- list()
    for(i in 1:mod$pdes){
      prior[[i]] <- list(dist="uniform")
    }
  }
  if(is.list(prior)){
    #HANDLE PRIORS
    for(i in 1:(mod$pdes)){
      # 1. Scale prior truncation bounds to BASS scale
      if(is.null(prior[[i]]$trunc)){
        prior[[i]]$trunc <- mod$range.des[,i]
      }
      prior[[i]]$trunc <- scale_range(prior[[i]]$trunc, mod$range.des[,i])

      # 2. Handle each distribution type separately
      distribution = prior[[i]]$dist
      if(distribution == "normal"){
        if(is.null(prior[[i]]$weights)){
          num_mix <- length(prior[[i]]$mean)
          prior[[i]]$weights <- rep(1/num_mix, num_mix)
        }
        prior[[i]]$mean <- scale_range(prior[[i]]$mean, mod$range.des[,i])
        prior[[i]]$sd <- prior[[i]]$sd/(mod$range.des[2,i] - mod$range.des[1,i])
        prior[[i]]$z <- pnorm((prior[[i]]$trunc[2]-prior[[i]]$mean)/prior[[i]]$sd) - pnorm((prior[[i]]$trunc[1]-prior[[i]]$mean)/prior[[i]]$sd)
        cc <- sum(prior[[i]]$weights*prior[[i]]$z)
        prior[[i]]$weights <- prior[[i]]$weights/cc
        # DF: prior[[i]]$z
        # change weights with truncation
        # divide by cc instead to keep the same prior shape
        # does the truncation change the distribution shape in the non-truncated regions??

        #Check for extrapolation
        qq <- qnorm(c(0.001, 0.999), prior[[i]]$mean, prior[[i]]$sd)
        if((qq[1] < 0 - 0.1) | (qq[2] > 1 + 0.1)){
          warning('You are asking the emulator to extrapolate. This is not reccomended.')
        }
      }
      if(distribution == "uniform"){
        prior[[i]]$weights = 1
        #Check for extrapolation
        if(((prior[[i]]$trunc[1] < 0 - 0.1) | (prior[[i]]$trunc[2] > 1 + 0.1))){
          warning('You are asking the emulator to extrapolate. This is not reccomended.')
        }
      }
      if(distribution == "beta"){
        prior[[i]]$weights = 1
        if(abs(min(mod$xx.des[,i]) - 0) > 0.01 | abs(max(mod$xx.des[,i]) - 1) > 0.01){
          warning("For dist=\"beta\", data should be scaled to (0, 1) range BEFORE fitting bass model")
        }
      }
      if(distribution == "gamma"){
        prior[[i]]$weights = 1
        if(abs(min(mod$xx.des[,i]) - 0) > 0.01){
          warning("For dist=\"gamma\", data should be shifted to have a minimum of 0 BEFORE fitting bass model")
        }
        prior[[i]]$rate = prior[[i]]$rate * prior[[i]]$trunc[2]
      }
    }

    #browser()
    # Make prior into a function
    rprior <- get_prior_generator(prior)
  }else{
    rprior <- prior
  }

  # Start monte carlo loop
  X_mc <- t(sapply(1:n_mc, function(xx) rprior()))
  grads1 <- gradient_bass(mod, X_mc, mcmc.use=mcmc.use[,1])
  grads2 <- gradient_bass(mod2, X_mc, mcmc.use=mcmc.use[,2])

  n_mcmc <- dim(grads1)[1]
  p <- dim(grads1)[3]
  C_out <- list()
  A_tform <- diag(1/apply(mod$range.des, 2, diff))
  for(i in seq_len(n_mcmc)){
    C_tmp <- matrix(0, nrow=p, ncol=p)
    for(j in seq_len(n_mc)){
      C_tmp <- C_tmp + tcrossprod(grads1[i,j,], grads2[i,j,]) / n_mc
    }

    if(use_native_scale){
      # Transform back to native space
      C_tmp <- t(A_tform)%*%C_tmp%*%A_tform
    }
    C_out[[i]] <- C_tmp
  }
  if(n_mcmc == 1){
    return(C_out[[1]])
  }
  return(C_out)
}

