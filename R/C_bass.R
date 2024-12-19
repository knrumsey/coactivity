#' Estimate the Constantine Matrix with BASS
#'
#' Closed form estimator of the C matrix using a BASS model
#'
#' @param mod a fitted BASS model. The output of the bass() or bassPCA() functions.
#' @param prior a list, like one returned by the \code{build_prior()} function. See the documentation for details.
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
C_bass <- function(mod, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.use=NULL, func.true=NULL, pca_method=2, verbose=FALSE){
  if(class(mod) == "list"){
    nmods <- length(mod)
    res <- list()
    for(i in seq_len(nmods)){
      res[[i]] <- C_bass(mod[[i]], prior, mcmc.use, use_native_scale, func.use[i], pca_method, verbose)
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
          res[[i]] <- C_bass_univariate(mod_t[[i]], prior, mcmc.use, use_native_scale, verbose)
        }
      }else{
        res <- C_bass_univariate(mod_t, prior, mcmc.use, use_native_scale, verbose)
      }
      return(res)
    }else{
      res <- C_bass_univariate(mod, prior, mcmc.use, use_native_scale, verbose)
      return(res)
    }
  }else if("bassBasis" %in% class(mod)){
    if(pca_method == 1){
      C_bass_pca_v1(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose)
    }else{
      C_bass_pca_v2(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose)
    }

  }else{
    return("mod should be a fitted bass model")
  }
}


get_u <- function(t, U, k){
  return(approx(seq(0,1,length.out=nrow(U)),
         U[,k],
         t)$y)
}

C_bass_pca_v1 <- function(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
  return(TRUE)
}

C_bass_pca_v2 <- function(mod, prior, mcmc.use, func.use, func.true, use_native_scale, verbose){
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
    res[[cnt]] <- C_bass_univariate(mod_curr, prior, mcmc.use, use_native_scale, verbose)
    cnt <- cnt + 1
  }

  if(length(func.use) == 1){
    return(res[[1]])
  }else{
    return(res)
  }
}


C_bass_univariate <- function(mod, prior = NULL, mcmc.use=NULL, use_native_scale=FALSE, verbose=TRUE){
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

  # Make constantine matrix
  Cf_post <- list()
  # Get transformation matrix
  A_tform <- diag(1/apply(mod$range.des, 2, diff))
  for(r in seq_along(mcmc.use)){
    #Compute only the stuff we will need for every iteration
    rr <- mcmc.use[r]
    mod_number_new <- mod$model.lookup[rr]
    coeff      <- mod$beta[rr,]
    coeff      <- matrix(coeff[!is.na(coeff)][-1], nrow=1)
    M_new <- length(coeff)

    compute_flag <- FALSE
    if(r == 1){
      compute_flag <- TRUE
    }else{
      if(mod_number != mod_number_new){
        compute_flag <- TRUE
      }
    }

    if(compute_flag){
      mod_number <- mod_number_new
      M          <- M_new
      signs      <- mod$signs.des[mod_number, 1:M, ]
      indic      <- mod$vars.des[mod_number, 1:M, ]
      knots      <- mod$knotInd.des[mod_number, 1:M, ]

      if(M==0){
        Cf <- matrix(0, nrow=mod$pdes, ncol=mod$pdes)
      }else{
        if(M==1){
          signs <- matrix(signs, nrow=1, ncol=length(signs))
          indic <- matrix(indic, nrow=1, ncol=length(indic))
          knots <- matrix(knots, nrow=1, ncol=length(knots))
        }
        #if(gbass_flag){
        #  knots    <- mod$knots.des[mod_number, 1:M, ]
        #}else{
        #  knots    <- mod$knotInd.des[mod_number, 1:M, ]
        #}
        # Initalize arrays
        C <- A <- B <- I1 <- I2 <- I3 <- array(NA, dim=c(mod$pdes, M, M))

        for(i in 1:mod$pdes){
          prior_i <- prior[[i]]
          num_c <- ncol(signs)
          v <- apply(indic, 1, function(zz) match(i, zz))
          u <- !is.na(v)
          s <- apply(cbind(signs, v), 1, function(zz) zz[zz[num_c + 1]])
          if(gbass_flag || bassfunc_flag){
            t <- apply(cbind(knots, v), 1, function(zz) zz[zz[num_c + 1]])
          }else{
            t <- Xt[apply(cbind(knots, v), 1, function(zz) zz[zz[num_c + 1]]), i]
          }

          # If this comes from BASS (rather than GBASS with gm2bm)
          # then we need to account for Devin's g-scaling-factors
          #if(!gbass_flag){
          d <- 1/((s + 1)/2 - s*t)
          s <- s*d
          #}

          #Handle NA cases
          s[is.na(s)] <- 1
          t[is.na(t)] <- -Inf

          # Get integration constants
          c1 <- tcrossprod(s)
          C[i,,] <- c1

          a <- bp <- i1 <- i2 <- i3 <- matrix(0, M, M)
          for(m1 in 1:M){
            for(m2 in m1:M){
              um  <- u[c(m1, m2)]
              ssm <- s[c(m1, m2)]
              sm  <- sign(ssm)
              tm  <- t[c(m1, m2)]

              if(sm[1] == 1 & sm[2] == 1){
                a[m1,m2]  <- a[m2,m1]  <- max(tm)
                bp[m1,m2] <- bp[m2,m1] <- Inf
              }
              if(sm[1] == 1 & sm[2] == -1){
                a[m1,m2]  <- a[m2,m1]  <- tm[1]
                bp[m1,m2] <- bp[m2,m1] <- tm[2]
              }
              if(sm[1] == -1 & sm[2] == 1){
                a[m1,m2]  <- a[m2,m1]  <- tm[2]
                bp[m1,m2] <- bp[m2,m1] <- tm[1]
              }
              if(sm[1] == -1 & sm[2] == -1){
                a[m1,m2]  <- a[m2,m1]  <- -Inf
                bp[m1,m2] <- bp[m2,m1] <- min(tm)
              }
              aa <- a[m1,m2]
              bb <-  max(a[m1,m2], bp[m1,m2])

              # Compute integrals
              #E0 <- Efunc(0, aa, bb, prior_i)
              #E1 <- Efunc(1, aa, bb, prior_i)
              #E2 <- Efunc(2, aa, bb, prior_i)

              E0 <- XI_FUNC(0, aa, bb, prior_i)
              E1 <- XI_FUNC(1, aa, bb, prior_i)
              E2 <- NULL #Only compute when needed

              if(is.nan(E1) | is.nan(E0)){
                browser()
              }

              # Start with I3
              if(um[1] == 0 | um[2] == 0){
                i3[m1,m2] <- i3[m2,m1] <- 0
              }else{
                i3[m1,m2] <- i3[m2,m1] <- E0
              }

              #Next, I1 and I2
              if(um[1] == 0){
                if(um[2] == 0){
                  i1[m1,m2] <- 0
                  i1[m2,m1] <- 0
                  i2[m1,m2] <- i2[m2,m1] <- E0
                }else{
                  i1[m1,m2] <- 0
                  i1[m2,m1] <- E0
                  i2[m1,m2] <- i2[m2,m1] <- E1 - tm[2]*E0
                }
              }else{
                if(um[2] == 0){
                  i2[m1,m2] <- i2[m2,m1] <- E1 - tm[1]*E0
                  i1[m1,m2] <- E0
                  i1[m2,m1] <- 0
                }else{
                  E2 <- XI_FUNC(2, aa, bb, prior_i)
                  i1[m1,m2] <- E1 - tm[2]*E0
                  i1[m2,m1] <- E1 - tm[1]*E0
                  i2[m1,m2] <- i2[m2,m1] <- E2 - (tm[1] + tm[2])*E1 + tm[1]*tm[2]*E0
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
        }
        #Add signs to the I's
        I1 <- I1*C
        I2 <- I2*C
        I3 <- I3*C
      }
    }
    if(M > 0){
      #Reconstruct Constantine matrix
      Cf <- matrix(NA, nrow=mod$pdes, ncol=mod$pdes)
      for(i in 1:mod$pdes){
        for(j in 1:i){
          if(i == j){
            cij_curr <- crossprod(coeff)*I3[i,,]
            for(k in (1:mod$pdes)[-c(i)]){
              cij_curr <- cij_curr * I2[k,,]
            }
            Cf[i,j] <- sum(cij_curr)
          }else{
            cij_curr <- crossprod(coeff)*I1[i,,]*t(I1[j,,])
            for(k in (1:mod$pdes)[-c(i,j)]){
              cij_curr <- cij_curr * I2[k,,]
            }
            Cf[i,j] <- Cf[j,i] <- sum(cij_curr)
          }
          # Matrix-y way of doing things (kinda)
          # This might save a little time when p is large, but the way I chose to code it (for time savings) has divide by zero problems if not careful
          # if(mod$pdes > 30){
          #   I2_prod <- matrix(1, M, M)
          #   for(k in (1:mod$pdes)){
          #     I2_prod <- I2_prod * I2[k,,]
          #   }
          #   if(i == j){
          #     Cf[i,i] <- coeff%*%(I3[i,,] * I2_prod / (I2[i,,]+ 1e-12))%*%t(coeff)
          #   }else{
          #     Cf[i,j] <- Cf[j,i] <- coeff%*%(I1[i,,] * t(I1[j,,]) * I2_prod / (I2[i,,] + 1e-12) / (I2[j,,] + 1e-12))%*%t(coeff)
          #   }
          # }else{do it another way}
          #A third way to do it
          # I2_prod <- matrix(1, M, M)
          # if(i == j){
          #   for(k in (1:mod$pdes)[-i]){
          #     I2_prod <- I2_prod * I2[k,,]
          #   }
          #   Iii <- I3[i,,] * I2_prod
          #   Cf[i,i] <- coeff%*%Iii%*%t(coeff)
          # }else{
          #   for(k in (1:mod$pdes)[-c(i,j)]){
          #     I2_prod <- I2_prod * I2[k,,]
          #   }
          #   Iij <- I1[i,,] * t(I1[j,,]) * I2_prod
          #   Cf[i,j] <- Cf[j,i] <- coeff%*%Iij%*%t(coeff)
          # }
        }
      }
      if(use_native_scale){
        # Transform back to native space
        Cf <- t(A_tform)%*%Cf%*%A_tform
      }
    }
    if(verbose & length(mcmc.use > 5) & r %in% round(seq(0, length(mcmc.use), length.out=5))[-1]){
      cat("Processing MCMC iteration", mcmc.use[r], myTimestamp(), "\n")
    }
    Cf_post[[r]] <- Cf
  }
  class(Cf_post) <- "ConstantineMatrix"
  if(length(Cf_post) == 1){
    return(Cf_post[[1]])
  }
  return(Cf_post)
}




# Computes the truncated moments with respect to the prior
XI_FUNC <- function(pow, a, b, prior){
  astar <- max(a, prior$trunc[1])
  bstar <- max(astar, min(b, prior$trunc[2]))

  if(bstar > astar){
    # Number of components in mixture
    L <- length(prior$weights)
    res <- 0
    for(ell in 1:L){
      distribution <- prior$dist
      if(distribution == "uniform"){
        res <- res + XI_FUNC_UNIF(pow, astar, bstar, prior, ell)*prior$weights[ell]
      }
      if(distribution == "normal"){
        res <- res + XI_FUNC_TNORM_SVW(pow, astar, bstar, prior, ell)*prior$weights[ell]
      }
      if(distribution == "beta"){
        res <- res + XI_FUNC_BETA(pow, astar, bstar, prior, ell)*prior$weights[ell]
      }
      if(distribution == "gamma"){
        res <- res + XI_FUNC_GAMMA(pow, astar, bstar, prior, ell)*prior$weights[ell]
      }
    }
    return(res)
  }else{
    return(0)
  }
}

XI_FUNC_TNORM <- function(pow, a, b, prior, ell){
  mu    <- prior$mean[ell]
  sigma <- prior$sd[ell]
  tau0  <- prior$trunc[1]
  tau1  <- prior$trunc[2]
  AA    <- (a-mu)/sigma
  BB    <- (b-mu)/sigma
  if(is.infinite(AA)) AA <- sign(AA)*1e4*sigma
  if(is.infinite(BB)) BB <- sign(BB)*1e4*sigma
  Z0    <- (pnorm(BB) - pnorm(AA))/(pnorm((tau1 - mu)/sigma) - pnorm((tau0-mu)/sigma))
  term  <- 1
  if(pow == 1){
    Z1 <- (dnorm(AA) - dnorm(BB))/(pnorm(BB) - pnorm(AA))
    term <- mu + sigma*Z1
  }
  if(pow == 2){
    Z1 <- (dnorm(AA) - dnorm(BB))/(pnorm(BB) - pnorm(AA))
    Z2 <- 1 - (AA*dnorm(AA) - BB*dnorm(BB))/(pnorm(BB) - pnorm(AA))
    term <- mu^2 + 2*sigma*mu*Z1 + sigma^2*Z2
  }
  return(Z0*term)
}


XI_FUNC_TNORM_SVW <- function(pow, a, b, prior){
  mu    <- prior$mean[ell]
  sigma <- prior$sd[ell]
  tau0  <- prior$trunc[1]
  tau1  <- prior$trunc[2]
  AA    <- (a-mu)/sigma
  BB    <- (b-mu)/sigma
  T1    <- (tau1-mu)/sigma
  T0    <- (tau0-mu)/sigma
  if(is.infinite(AA)) AA <- sign(AA)*1e6*sigma
  if(is.infinite(BB)) BB <- sign(BB)*1e6*sigma

  Delta1 <- pnorm(BB) - pnorm(AA)
  Delta2 <- pnorm(T1) - pnorm(T0)

  term  <- 0
  if(pow == 1){
    Z1     <- dnorm(AA) - dnorm(BB)
    term   <- sigma*Z1
  }
  if(pow == 2){
    Z1     <- dnorm(AA) - dnorm(BB)
    Z2     <- Delta1 - AA*dnorm(AA) - BB*dnorm(BB)
    term   <- 2*sigma*mu*Z1 + sigma^2*Z2
  }
  res <- Delta1*mu^pow/Delta2 + term/Delta2
  return(res)
}


XI_FUNC_UNIF <- function(pow, a, b, prior, ell){
  tau0  <- prior$trunc[1]
  tau1  <- prior$trunc[2]
  res <- (b^(pow+1) - a^(pow+1))/((tau1-tau0)*(pow+1))
  return(res)
}

XI_FUNC_BETA <- function(pow, a, b, prior, ell){
  shape1 <- prior$shape1[ell]
  shape2 <- prior$shape2[ell]
  num    <- zipfR::Ibeta(b, shape1 + pow, shape2) - zipfR::Ibeta(a, shape1 + pow, shape2)
  den    <- beta(shape1, shape2)
  return(num/den)
}

XI_FUNC_GAMMA <- function(pow, a, b, prior, ell){
  alpha <- prior$scale[ell]
  beta  <- prior$rate[ell]
  num   <- zipfR::Igamma(alpha+pow, beta*b) - zipfR::Igamma(alpha+pow, beta*a)
  den   <- beta^pow*gamma(alpha)
  return(num/den)
}











