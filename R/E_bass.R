#' Expected value of a BASS model
#'
#' Closed form estimator of the expected value of a BASS function
#'
#' @param mod a fitted BASS model. The output of the bass() or bassPCA() functions.
#' @param prior a list, like one returned by the \code{build_prior()} function. See the documentation for details.
#' @param mcmc.use a vector of indices telling which mcmc draws to use
#' @param func.use a vector indicating which values of the functional variable to compute C for, if applicable
#' @param verbose Doesn't do anything currently.
#' @return A list representing the posterior distribution of the Constantine matrix.
#' @details Returns the expected value of a BASS model.
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
#' mod <- bass(XX, y1)
#' E <- E_bass(mod)
#' @export
E_bass <- function(mod, prior=NULL, mcmc.use=NULL, func.use=NULL, verbose=FALSE){
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
          res[[i]] <- E_bass_univariate(mod_t[[i]], prior, mcmc.use)
        }
      }else{
        res <- E_bass_univariate(mod_t, prior, mcmc.use)
      }
      return(res)
    }else{
      res <- E_bass_univariate(mod, prior, mcmc.use)
      return(res)
    }
  }else if("bassBasis" %in% class(mod)){
    if(is.null(func.use)){
      func.use <- seq(0, 1, length.out=nrow(mod$dat$xx))
    }
    # Convert to univariate bass models
    mod_t <- bassPCA2bass_fixed_t(mod, func.use, verbose)
    if(length(func.use) > 1){
      res <- list()
      for(i in seq_along(func.use)){
        res[[i]] <- E_bass_univariate(mod_t[[i]], prior, mcmc.use)
      }
    }else{
      res <- E_bass_univariate(mod_t, prior, mcmc.use)
    }
    return(res)
  }else{
    return("mod should be a fitted bass model")
  }
}




E_bass_univariate <- function(mod, prior = NULL, mcmc.use=NULL){
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

  # Make Z vector
  Ef_post <- list()
  # Get transformation vector
  A_tform <- 1/apply(mod$range.des, 2, diff)
  for(r in 1:length(mcmc.use)){
    #Compute only the stuff we will need for every iteration
    rr <- mcmc.use[r]
    mod_number_new <- mod$model.lookup[rr]
    coeff      <- mod$beta[rr,]
    intercept  <- coeff[1]
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
      if(M==1){
        signs <- matrix(signs, nrow=1, ncol=length(signs))
        indic <- matrix(indic, nrow=1, ncol=length(indic))
        knots <- matrix(knots, nrow=1, ncol=length(knots))
      }

      # Initalize arrays
      A <- B <- I5 <- Eim <- array(NA, dim=c(mod$pdes, M))
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
        # NOTE: Need to comment out, i don't remember why
        d <- 1/((s + 1)/2 - s*t)
        s <- s*d
        #}

        #Handle NA cases
        s[is.na(s)] <- 1
        t[is.na(t)] <- -Inf

        a <- b <- i5 <- matrix(0, M)
        #browser()
        for(m in 1:M){
          um <- u[m]
          ssm <- s[m]
          sm <- sign(ssm)
          tm <- t[m]

          a[m] <- ifelse(sm==1, tm, -Inf)
          b[m] <- ifelse(sm==1, Inf, tm)

          # Compute truncated moments
          E0  <- XI_FUNC(0, a[m], b[m], prior_i)
          E1  <- XI_FUNC(1, a[m], b[m], prior_i)

          if(is.nan(E1) | is.nan(E0)){
            browser()
          }
          # Compute integrals
          i5[m] <- ifelse(um == 0, 1, ssm*(E1 - tm*E0))
        }

        A[i,]  <- a
        B[i,]  <- b
        I5[i,] <- i5
      }
    }
    #Reconstruct Constantine matrix
    zi_curr <- coeff
    for(k in 1:mod$pdes){
      zi_curr <- zi_curr * I5[k,]
    }
    Ef <- sum(zi_curr) + intercept

    Ef_post[[r]] <- Ef
  }
  if(length(Ef_post) == 1){
    return(Ef_post[[1]])
  }
  return(Ef_post)
}
