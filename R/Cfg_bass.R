#' Estimate the Constantine Matrix with BASS
#'
#' Closed form estimator of the C matrix using a BASS model
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
#' C1_fg <- Cfg_bass(mod1_f, mod1_g)
#'
#' # ===========================================
#' #      CASE 2: Augmented BASS (fixed t)
#' # ===========================================
#' mod2_full_f <- bass(XX, yfunc_f, xx.func = xfunc)
#' mod2_full_g <- bass(XX, yfunc_g, xx.func = xfunc)
#' mod2_f <- bassfunc2bass_fixed_t(mod2_full_f, 0.5)
#' mod2_g <- bassfunc2bass_fixed_t(mod2_full_g, 0.5)
#' C2_fg <- Cfg_bass(mod2_f, mod2_g)
#' C2b_fg <- Cfg_bass(mod2_full_f, mod2_full_g, func.use = 0.5)
#'
#' # ===========================================
#' #      CASE 3: PCA BASS (fixed t)
#' # ===========================================
#' mod3_full_f <- bassPCA(XX, yfunc_f)
#' mod3_full_g <- bassPCA(XX, yfunc_g)
#' mod3_f <- bassPCA2bass_fixed_t(mod3_full_f, 0.5)
#' mod3_g <- bassPCA2bass_fixed_t(mod3_full_g, 0.5)
#' C3_fg <- Cfg_bass(mod3
#' @export
Cfg_bass <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.use=NULL, func.true=NULL, pca_method=2, verbose=FALSE){
  # HANDLE CASE WHERE mod1 or mod2 IS A LIST
  if(class(mod1) == "list" && class(mod2) == "list"){
    nmods <- length(mod1)
    res <- list()
    for(i in seq_len(nmods)){
      res[[i]] <- Cfg_bass(mod1[[i]], mod2[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
    }
    return(res)
  }else{
    if(class(mod1) == "list"){
      nmods <- length(mod1)
      res <- list()
      for(i in seq_len(nmods)){
        res[[i]] <- Cfg_bass(mod1[[i]], mod2, prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
      }
      return(res)
    }
    if(class(mod1) == "list"){
      nmods <- length(mod2)
      res <- list()
      for(i in seq_len(nmods)){
        res[[i]] <- Cfg_bass(mod1, mod2[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
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
    return(Cfg_bass_univariate(mod1, mod2, prior, mcmc.use, use_native_scale, verbose))
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
      res[[cnt]] <- Cfg_bass_univariate(mod1_t, mod2_t, prior, mcmc.use, use_native_scale, verbose)
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
        res[[cnt]] <- Cfg_bass_univariate(scal_mod, func_mod_t, prior, mcmc.use, use_native_scale, verbose)
      }else{
        res[[cnt]] <- Cfg_bass_univariate(func_mod_t, scal_mod, prior, mcmc.use, use_native_scale, verbose)
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
      return(Cfg_bass_pca_v1(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose))
    }else{
      return(Cfg_bass_pca_v2(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose))
    }
  }
}

Cfg_bass_pca_v1 <- function(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
  return(TRUE)
}

Cfg_bass_pca_v2 <- function(mod1, mod2, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
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

    res[[cnt]] <- Cfg_bass_univariate(curr1, curr2, prior, mcmc.use, use_native_scale, verbose)
    cnt <- cnt + 1
  }

  if(nfunc == 1){
    return(res[[1]])
  }else{
    return(res)
  }
}

Cfg_bass_univariate <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, verbose=FALSE){
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
      prior[[i]]$weights <- prior[[i]]$weights/cc # DF: prior[[i]]$z # change weights with truncation # divide by cc instead to keep the same prior shape# does the truncation change the distribution shape in the non-truncated regions??

      #Check for extrapolation
      qq <- qnorm(c(0.0005, 0.9995), prior[[i]]$mean, prior[[i]]$sd)
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

  # Make constantine matrix
  Cfg_post <- list()
  A_tform <- diag(1/apply(mod$range.des, 2, diff))
  Xt <- mod$xx.des
  for(r in 1:nrow(mcmc.use)){
    #Compute only the stuff we will need for every iteration
    rr1 <- mcmc.use[r,1]
    rr2 <- mcmc.use[r,2]

    mod_number_new <- mod$model.lookup[rr1]
    coeff          <- mod$beta[rr1,]
    coeff          <- matrix(coeff[!is.na(coeff)][-1], nrow=1)
    M_new          <- length(coeff)

    mod_number_new2 <- mod2$model.lookup[rr2]
    coeff2          <- mod2$beta[rr2,]
    coeff2          <- matrix(coeff2[!is.na(coeff2)][-1], nrow=1)
    M_new2          <- length(coeff2)

    compute_flag <- FALSE
    if(r == 1){
      compute_flag <- TRUE
    }else{
      if((mod_number != mod_number_new) | (mod_number2 != mod_number_new2)){
        compute_flag <- TRUE
      }
    }

    if(compute_flag){
      mod_number <- mod_number_new
      mod_number2 <- mod_number_new2
      M          <- M_new
      M2          <- M_new2

      if(any(c(M, M2) == 0)){
        Cfg <- matrix(0, nrow=mod$pdes, ncol=mod$pdes)
      }else{
        knots <- mod$knotInd.des[mod_number, 1:M, ]
        signs <- mod$signs.des[mod_number, 1:M, ]
        indic <- mod$vars.des[mod_number, 1:M, ]
        if(M==1){
          signs <- matrix(signs, nrow=1, ncol=length(signs))
          indic <- matrix(indic, nrow=1, ncol=length(indic))
          knots <- matrix(knots, nrow=1, ncol=length(knots))
        }
        knots2 <- mod2$knotInd.des[mod_number2, 1:M2, ]
        signs2 <- mod2$signs.des[mod_number2, 1:M2, ]
        indic2 <- mod2$vars.des[mod_number2, 1:M2, ]
        if(M2==1){
          signs2 <- matrix(signs2, nrow=1, ncol=length(signs2))
          indic2 <- matrix(indic2, nrow=1, ncol=length(indic2))
          knots2 <- matrix(knots2, nrow=1, ncol=length(knots2))
        }

        # Initalize arrays
        C <- A <- B <- I1 <- I2 <- I3 <- array(NA, dim=c(mod$pdes, M, M2))
        I1b <- array(NA, dim=c(mod$pdes, M2, M))
        for(i in 1:mod$pdes){
          prior_i <- prior[[i]]
          num_c <- ncol(signs)
          v <- apply(indic, 1, function(zz) match(i, zz))
          u <- !is.na(v)
          s <- apply(cbind(signs, v), 1, function(zz) zz[zz[num_c + 1]])
          if(gbass_flag1 || bassfunc_flag1){
            t <- apply(cbind(knots, v), 1, function(zz) zz[zz[num_c + 1]])
          }else{
            t <- Xt[apply(cbind(knots, v), 1, function(zz) zz[zz[num_c + 1]]), i]
          }

          num_c2 <- ncol(signs2)
          v2 <- apply(indic2, 1, function(zz) match(i, zz))
          u2 <- !is.na(v2)
          s2 <- apply(cbind(signs2, v2), 1, function(zz) zz[zz[num_c2 + 1]])
          if(gbass_flag2 || bassfunc_flag2){
            t2 <- apply(cbind(knots2, v2), 1, function(zz) zz[zz[num_c2 + 1]])
          }else{
            t2 <- Xt[apply(cbind(knots2, v2), 1, function(zz) zz[zz[num_c2 + 1]]), i]
          }
          # If this comes from BASS (rather than GBASS with gm2bm)
          # then we need to account for Devin's g-scaling-factors
          if(!gbass_flag1){
            d <- 1/((s + 1)/2 - s*t)
            s <- s*d
          }
          if(!gbass_flag2){
            d <- 1/((s2 + 1)/2 - s2*t2)
            s2 <- s2*d
          }
          #NOTE!!! DOES THE ABOVE WORK? CAN I HAVE A GBASS AND A BASS MODEL?

          #Handle NA cases
          s[is.na(s)] <- 1
          t[is.na(t)] <- -Inf

          s2[is.na(s2)] <- 1
          t2[is.na(t2)] <- -Inf

          # Get integration constants
          cc <- tcrossprod(s, s2)
          C[i,,] <- cc

          a <- bp <- i1 <- i2 <- i3 <- matrix(NA, M, M2)
          i1b <- a2 <- bp2 <- matrix(NA, M2, M)
          for(m1 in 1:M){
            for(m2 in 1:M2){
              #Main change for Cfg here (three lines)
              um <- c(u[m1], u2[m2])
              sm <- c(s[m1], s2[m2])
              tm <- c(t[m1], t2[m2])

              signm = sign(sm)

              if(signm[1] == 1 & signm[2] == 1){
                a[m1,m2]   <- max(tm)
                bp[m1,m2]  <- Inf
              }
              if(signm[1] == 1 & signm[2] == -1){
                a[m1,m2]   <- tm[1]
                bp[m1,m2]  <- tm[2]
              }
              if(signm[1] == -1 & signm[2] == 1){
                a[m1,m2]   <- tm[2]
                bp[m1,m2]  <- tm[1]
              }
              if(signm[1] == -1 & signm[2] == -1){
                a[m1,m2]  <- -Inf
                bp[m1,m2] <- min(tm)
              }
              aa <- a[m1,m2]
              bb <-  max(aa, bp[m1,m2])

              # Compute integrals
              E0 <- XI_FUNC(0, aa, bb, prior_i)
              E1 <- XI_FUNC(1, aa, bb, prior_i)
              E2 <- XI_FUNC(2, aa, bb, prior_i)

              # Start with I3
              if(um[1] == 0 | um[2] == 0){
                i3[m1,m2] <- 0
              }else{
                i3[m1,m2] <- E0
              }

              #Next, I1 and I2
              if(um[1] == 0){
                if(um[2] == 0){
                  i1[m1,m2]  <- 0
                  i1b[m2,m1] <- 0
                  i2[m1,m2]  <- E0
                }else{
                  i1[m1,m2]  <- 0
                  i1b[m2,m1] <- E0
                  i2[m1,m2]  <- E1 - tm[2]*E0
                }
              }else{
                if(um[2] == 0){
                  i1[m1,m2]  <- E0
                  i1b[m2,m1] <- 0
                  i2[m1,m2]  <- E1 - tm[1]*E0
                }else{
                  i1[m1,m2]  <- E1 - tm[2]*E0
                  i1b[m2,m1] <- E1 - tm[1]*E0
                  i2[m1,m2]  <- E2 - (tm[1] + tm[2])*E1 + tm[1]*tm[2]*E0
                }
              }
            }
          }
          b <- pmax(a, bp)
          A[i,,] <- a
          B[i,,] <- b
          I1[i,,] <- i1
          I2[i,,] <- i2
          I3[i,,] <- i3
          I1b[i,,] <- i1b
        }
        #Add signs to the I's
        I1 <- I1*C
        I2 <- I2*C
        I3 <- I3*C
        for(iii in 1:mod$pdes){
          I1b[iii,,] <- I1b[iii,,] * t(C[iii,,])
        }
      }
    } # Done computing

    if(all(c(M, M2) > 0)){
      #Reconstruct Constantine matrix
      Cfg <- matrix(NA, nrow=mod$pdes, ncol=mod$pdes)
      for(i in 1:mod$pdes){
        for(j in 1:mod$pdes){
          # #Naive way for now
          cij_curr <- 0
          if(i != j){
            for(m1 in 1:M){
              for(m2 in 1:M2){
                term <- coeff[m1]*coeff2[m2]*I1[i,m1,m2]*I1b[j,m2,m1]
                for(k in (1:mod$pdes)[-c(i,j)]){
                  term <- term * I2[k,m1,m2]
                }
                cij_curr <- cij_curr + term
              }
            }
          }else{
            for(m1 in 1:M){
              for(m2 in 1:M2){
                term <- coeff[m1]*coeff2[m2]*I3[i,m1,m2]
                for(k in (1:mod$pdes)[-i]){
                  term <- term * I2[k,m1,m2]
                }
                cij_curr <- cij_curr + term
              }
            }
          }
          Cfg[i,j] <- cij_curr
        }
      }
      if(use_native_scale == FALSE){
        # Transform back to native space
        Cfg <- t(A_tform)%*%Cfg%*%A_tform
      }
    } #End non-trivial Cfg calculation
    if(verbose & nrow(mcmc.use) > 5 & r %in% round(seq(0, nrow(mcmc.use), length.out=5))[-1]){
      cat("Processing MCMC iteration (", mcmc.use[r,1], ", ", mcmc.use[r,2], ") ", myTimestamp(), "\n", sep="")
    }
    #Cfg_post[[r]] <- Cfg
    Cfg_post[[paste0(rr1, ", ", rr2)]] <- Cfg
  }
  class(Cfg_post) <- "CoConstantineMatrix"
  if(length(Cfg_post) == 1){
    return(Cfg_post[[1]])
  }
  return(Cfg_post)
}
