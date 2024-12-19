#' Estimate the Augmented Constantine Matrix
#'
#' Closed form estimator of the C_aug matrix using a functional BASS model
#'
#' @param mod a fitted (functional) BASS model. The output of the bass() or bassPCA() functions.
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
C_aug_bass <- function(mod, prior=NULL, mcmc.use=NULL, use_native_scale=FALSE, func.true=NULL, verbose=FALSE){
  if("bass" %in% class(mod)){
    if(mod$pfunc == 0){
      stop("mod must be a functional BASS model.")
    }
    mod_new <- bassfunc2bass(mod)
    res <- C_bass(mod_new, prior, mcmc.use, use_native_scale, verbose=verbose)
  }else if(class(mod) == "bassBasis"){
    nfunc <- nrow(mod$dat$basis)
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
    res <- C_aug_bassPCA(mod, prior, rho_t, mcmc.use, use_native_scale, func.true, verbose)
  }else{
    stop("mod must be a functional BASS model")
  }
  return(res)
}

C_aug_bassPCA <- function(mod, prior, rho_t, mcmc.use, use_native_scale, func.true, verbose){
  if(is.null(mcmc.use)){
    mcmc.use <- length(mod$mod.list[[1]]$nbasis)
  }
  nmc <- length(mcmc.use)
  U <- mod$dat$basis
  mod_list <- mod$mod.list
  K <- length(mod_list)
  p <- mod_list[[1]]$pdes
  if(nmc == 1){
    Caug <- matrix(0, nrow=p+1, ncol=p+1)
  }else{
    Caug <- list()
    for(i in seq_len(nmc)){
      Caug[[i]] <- matrix(0, nrow=p+1, ncol=p+1)
    }
  }
  for(k1 in 1:K){
    for(k2 in 1:K){
      obj <- Cfg_aug_bass_univariate(mod_list[[k1]], mod_list[[k2]], prior, mcmc.use, verbose)

      u1 <- U[,k1]
      u2 <- U[,k2]
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

Cfg_aug_bass_univariate <- function(mod1, mod2, prior=NULL, mcmc.use=NULL, verbose=FALSE){
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
  Cfg_post <- Lfg_post <- Lgf_post <- lfg_post <- list()
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

      # Get Lamba Vectors
      Lambda_fg <- Lambda_gf <- rep(NA, mod$pdes)
      for(j in 1:mod$pdes){
        curr_j_fg <- curr_j_gf <- 0
        if(j == 1){
          lambda_fg <- 0
        }
        for(m1 in 1:M){
          for(m2 in 1:M2){
            term_fg <- coeff[m1]*coeff2[m2]*I1b[j,m2,m1]
            term_gf <- coeff[m1]*coeff2[m2]*I1[j,m1,m2]
            if(j==1){
              term_0 <- coeff[m1]*coeff2[m2]*I2[j,m1,m2]
            }
            for(k in (1:mod$pdes)[-j]){
              term_fg <- term_fg * I2[k,m1,m2]
              term_gf <- term_gf * I2[k,m1,m2]
              if(j == 1){
                term_0 <- term_0 * I2[k,m1,m2]
              }
            }
            curr_j_fg <- curr_j_fg + term_fg
            curr_j_gf <- curr_j_gf + term_gf
            if(j == 1){
              lambda_fg <- lambda_fg + term_0
            }
          }
        }
        Lambda_fg[j] <- curr_j_fg
        Lambda_gf[j] <- curr_j_gf
      }
      # DONT DO RESCALING HERE ANYMORE

    } #End non-trivial Cfg calculation
    if(verbose & nrow(mcmc.use) > 5 & r %in% round(seq(0, nrow(mcmc.use), length.out=5))[-1]){
      cat("Processing MCMC iteration (", mcmc.use[r,1], ", ", mcmc.use[r,2], ") ", myTimestamp(), "\n", sep="")
    }
    #Cfg_post[[r]] <- Cfg
    Cfg_post[[paste0(rr1, ", ", rr2)]] <- Cfg
    Lfg_post[[paste0(rr1, ", ", rr2)]] <- Lambda_fg
    Lgf_post[[paste0(rr1, ", ", rr2)]] <- Lambda_gf
    lfg_post[[paste0(rr1, ", ", rr2)]] <- lambda_fg
  }
  if(length(Cfg_post) == 1){
    out <- list(Cfg=Cfg_post[[1]], Lfg=Lfg_post[[1]], Lgf=Lgf_post[[1]], lfg=lfg_post[[1]])
    return(out)
  }
  out <- list(Cfg=Cfg_post, Lfg=Lfg_post, Lgf=Lgf_post, lfg=lfg_post)
  return(out)
}


compute_derivative <- function(u, t) {
  n <- length(u)
  u_deriv <- numeric(n)

  # One-sided differences for endpoints
  u_deriv[1] <- (u[2] - u[1]) / (t[2] - t[1])
  u_deriv[n] <- (u[n] - u[n-1]) / (t[n] - t[n-1])

  # Central differences for interior points
  for (i in 2:(n-1)) {
    u_deriv[i] <- (u[i+1] - u[i-1]) / (t[i+1] - t[i-1])
  }

  return(u_deriv)
}

# Function to compute integrals using the trapezoidal rule
trapezoidal <- function(t, f) {
  sum((t[-1] - t[-length(t)]) * (f[-1] + f[-length(f)]) / 2)
}


# Scalar multiple of a * every component of lst
# If lst is not a list object, then return a * lst
list_multiply <- function(a, lst){
  if("list" %in% class(lst)){
    res <- list()
    for(i in seq_along(lst)){
      res[[i]] <- a * lst[[i]]
    }
    return(res)
  }else{
    return(a * lst)
  }
}

