#' Extract scalar BASS model from a functional BASS model for fixed functional variable values
#'
#' Extracts the BASS model(s) corresponding to fixed values of the functional variable.
#' For each specified \( t \)-value, calculates the contribution of the functional variable
#' to the tensor product basis functions and returns a modified BASS model.
#'
#' @param bassfunc A BASS model object (with functional variables).
#' @param func.use A numeric vector of fixed values for the functional variable \( t \).
#' @param verbose logical; should progress be displayed?
#' @return A list of BASS models for each \( t \)-value (or a single model if length(func.use) == 1).
#' @examples
#' # Simulate bass models for each of the three cases
#' f<-function(x){
#'   10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
#' }
#' n<-500 # number of observations
#' nfunc<-50 # size of functional variable grid
#' xfunc<-seq(0,1,length.out=nfunc) # functional grid
#' x<-matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
#' X<-cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
#' y<-matrix(f(X),nrow=n)
#'
#' # Scalar response
#' X <- lhs::randomLHS(500, 5)
#' y <- f(X)
#' bassfunc <- bass(x, y, xx.func=xfunc)
#'
#' # Extract the BASS model for t = 0.5
#' mod_fixed <- funcbass2bass_fixed_t(mod_func, func.use = 0.5)
#'
#' # Extract models for multiple t values
#' mods_fixed <- funcbass2bass_fixed_t(mod_func, func.use = c(0.5, 0.75))
#' @export
bassfunc2bass_fixed_t <- function(bassfunc, func.use, verbose=FALSE) {
  if (bassfunc$pfunc == 0) {
    stop("The provided BASS model does not include a functional variable.")
  }

  # Initialize results
  results <- list()

  # Redefine functional variable numbers
  bassfunc$vars.func <- bassfunc$vars.func + bassfunc$pdes

  # Loop over each t value in func.use
  if (verbose & length(func.use) > 1) {
    nfunc <- length(func.use)
    verbose_vec <- round(seq(1, nfunc, length.out=5))
    cat(nfunc, "models to convert", myTimestamp(), "\n")
  }
  cnt <- 1
  for(t_val in func.use){
    if(verbose & length(func.use) > 1){
      if (cnt %in% verbose_vec) {
        cat("Converting model", cnt, myTimestamp(), "\n")
      }
    }
    mod <- bassfunc # For returning later
    nm <- bassfunc$n.models
    maxBasis <- max(bassfunc$nbasis)

    # Copy data for manipulation
    n.int <- bassfunc$n.int.des + bassfunc$n.int.func
    signs <- abind_custom(bassfunc$signs.des, bassfunc$signs.func, along=3)
    vars <- abind_custom(bassfunc$vars.des, bassfunc$vars.func, along=3)
    beta <- bassfunc$beta[,1:(1+maxBasis)]

    # Extract knots
    # Initialize knots.des and knots.func arrays
    knots.des <- array(NA, dim = dim(bassfunc$knotInd.des))
    knots.func <- array(NA, dim = dim(bassfunc$knotInd.func))

    nint.des <- dim(knots.des)[3]
    nint.func <- dim(knots.func)[3]
    maxInt <- max(nint.des, nint.func)
    for(i in seq_len(dim(bassfunc$knotInd.des)[1])){  # Loop over models
      for(j in seq_len(dim(bassfunc$knotInd.des)[2])){  # Loop over basis functions
        for(k in seq_len(maxInt)){  # Loop over interaction terms
          # For design variables
          if(k <= nint.des){
            ii.des <- bassfunc$knotInd.des[i, j, k]
            jj.des <- bassfunc$vars.des[i, j, k]
            if(!is.na(ii.des) && !is.na(jj.des)){
              knots.des[i, j, k] <- bassfunc$xx.des[ii.des, jj.des]
            }
          }

          # For functional variables
          if(k <= nint.func){
            ii.func <- bassfunc$knotInd.func[i, j, k]
            jj.func <- bassfunc$vars.func[i, j, k] - mod$pdes
            if(!is.na(ii.func) && !is.na(jj.func)){
              knots.func[i, j, k] <- bassfunc$xx.func[ii.func, jj.func]
            }
          }
        }
      }
    }
    knots <- abind_custom(knots.des, knots.func, along=3)
    knotInds <- abind_custom(bassfunc$knotInd.des, bassfunc$knotInd.func, along=3) # We want these for finding duplicates

    for(i in 1:nm){
      mcmc_iters <- which(mod$model.lookup == i)
      nbasis_curr <- mod$nbasis[mcmc_iters[1]]
      for(j in 1:nbasis_curr){
        vars_curr <- vars[i, j, ]
        signs_curr <- signs[i, j, ]
        knots_curr <- knots[i, j, ]
        knotInds_curr <- knotInds[i, j, ]

        # Check if functional variable is in the current basis function
        func_idx <- which(vars_curr == (bassfunc$pdes + bassfunc$pfunc))

        if(length(func_idx) > 0){
          # Evaluate the functional part at t_val
          func_sign <- signs_curr[func_idx]
          func_knot <- knots_curr[func_idx]

          # Devin's rescaling thing (I'm pretty sure we need this)
          d_scale <- 1/((func_sign + 1)/2 - func_sign*func_knot)
          func_sign <- d_scale * func_sign

          # Compute functional piece of the tensor product
          hij <- func_sign * (t_val - func_knot)
          chi <- as.numeric(hij > 0)

          # Try to make the generalize to more than one functional variable
          if(length(hij) > 1){
            hij <- prod(hij)
            chi <- prod(chi)
          }

          # Update the coefficient by including the functional contribution
          beta[mcmc_iters, j + 1] <- beta[mcmc_iters, j + 1] * chi  * hij

          # Remove the functional variable from the basis function
          #n.int[i,j] <- n.int[i,j] - 1
          vars[i, j, func_idx] <- NA
          signs[i, j, func_idx] <- NA
          knots[i, j, func_idx] <- NA
          knotInds[i, j, func_idx] <- NA

          # Shift NA values to the end
          # This code does nothing?
          # non_na_indices <- which(!is.na(vars[i, j, ]))
          # vars[i, j, seq_along(non_na_indices)] <- vars[i, j, non_na_indices]
          # vars[i, j, (length(non_na_indices) + 1):length(vars[i, j, ])] <- NA
          #
          # signs[i, j, seq_along(non_na_indices)] <- signs[i, j, non_na_indices]
          # signs[i, j, (length(non_na_indices) + 1):length(signs[i, j, ])] <- NA
          #
          # knots[i, j, seq_along(non_na_indices)] <- knots[i, j, non_na_indices]
          # knots[i, j, (length(non_na_indices) + 1):length(knots[i, j, ])] <- NA
          #
          # knotInds[i, j, seq_along(non_na_indices)] <- knotInds[i, j, non_na_indices]
          # knotInds[i, j, (length(non_na_indices) + 1):length(knotInds[i, j, ])] <- NA
          #
          # foo1 <- sum(is.na(signs[i,j,]))
          # foo2 <- sum(is.na(knots[i,j,]))
          # if(foo1 != foo2){
          #   print(paste0(i, ", ", j))
          # }
        }

        # Set n.int
        n.int[i,j] <- sum(!is.na(vars[i,j,]))
      }

      # Combine duplicate basis functions
      # After removing t, this is possible again
      model_matrix <- cbind(knotInds[i,1:nbasis_curr,], signs[i,1:nbasis_curr,], vars[i,1:nbasis_curr,])
      duplicates <- find_duplicate_groups(model_matrix)

      for (dup_group in duplicates) {
        # Check for intercept
        is_intercept <- all(is.na(vars[i,dup_group[1],]))
        if(is_intercept){
          # Check this is there are problems
          beta_mat <- matrix(beta[mcmc_iters, 1 + dup_group], nrow=length(mcmc_iters))
          beta[mcmc_iters, 1] <- beta[mcmc_iters, 1] + rowSums(beta_mat)
          to_remove <- dup_group
        }else{
          beta_mat <- matrix(beta[mcmc_iters, 1 + dup_group[-1]], nrow=length(mcmc_iters))
          beta[mcmc_iters, 1 + dup_group[1]] <- beta[mcmc_iters, 1 + dup_group[1]] + rowSums(beta_mat)
          to_remove <- dup_group[-1]
        }

        # Remove duplicate rows from model arrays
        # Mark rows as NA instead of removing
        knots[i, to_remove, ] <- NA
        knotInds[i, to_remove, ] <- NA
        signs[i, to_remove, ] <- NA
        vars[i, to_remove, ] <- NA
        n.int[i,to_remove] <- NA
        beta[mcmc_iters, 1 + to_remove] <- NA
      }

      # Get which columns to keep
      ikeep <- which(!is.na(knots[i,,1]))
      nkeep <- length(ikeep)
      inew <- seq_len(nkeep)

      # Shift NA columns to the end
      knots[i,inew,] <- knots[i,ikeep,]
      knots[i,-inew,] <- NA

      knotInds[i,inew,] <- knotInds[i,ikeep,]
      knotInds[i,-inew,] <- NA

      signs[i,inew,] <- signs[i,ikeep,]
      signs[i,-inew,] <- NA

      vars[i,inew,] <- vars[i,ikeep,]
      vars[i,-inew,] <- NA

      n.int[i,inew] <- n.int[i,ikeep]
      n.int[i,-inew] <- NA

      ikeep_beta <- c(1, 1 + ikeep)
      inew_beta  <- seq_len(nkeep + 1)
      beta[mcmc_iters, inew_beta] <- beta[mcmc_iters, ikeep_beta]
      beta[mcmc_iters, -inew_beta] <- NA

      mod$nbasis[mcmc_iters] <- nkeep
    }

    # Trim NA's off of beta
    maxBasis <- max(mod$nbasis) + 1
    if(any(!is.na(beta[,maxBasis+1]))){
      warning("Removing betas that we shouldnt?")
    }
    beta <- beta[,1:maxBasis]

    # Add new fields
    mod$beta <- beta
    mod$n.int.des <- n.int
    mod$signs.des <- signs
    mod$vars.des  <- vars
    mod$knots.des <- knots
    mod$knotInd.des <- knotInds
    mod$pdes <- bassfunc$pdes
    mod$range.des <- bassfunc$range.des
    mod$maxInt.des <- bassfunc$maxInt.des
    mod$cx <- bassfunc$cx
    if(mod$cat){
      mod$type <- "_des_cat"
    }else{
      mod$type <- "_des"
    }

    # Remove invalid fields
    mod$des.basis <- NULL
    mod$func.basis <- NULL
    mod$curr.list <- NULL

    # Remove functional fields
    mod$pfunc     <- 0
    mod$n.int.func <- NULL
    mod$signs.func <- NULL
    mod$vars.func <- NULL
    mod$knotInd.func <- NULL
    mod$range.func <- NULL
    mod$maxInt.func <- NULL
    mod$func      <- FALSE
    mod$wasfunc   <- TRUE   # Use this to flag (in C_bass) that knots are not indices.
    class(mod)    <- c("bass")

    # Update y, yhat, and yhat.mean
    if(bassfunc$pfunc == 1){
      mod$y_full <- mod$y
      y_new <- sapply(seq_len(nrow(mod$y)), function(ii) approx(mod$xx.func, mod$y[ii,], t_val))
      y_new <- unlist(y_new[2,])
      mod$y <- y_new

      if(!is.null(mod$yhat)){
        mod$yhat_full <- mod$yhat
        yhat_new <- array(NA, dim = c(dim(bassfunc$yhat)[1], dim(bassfunc$yhat)[2], length(func.use)))
        for(m in seq_len(dim(bassfunc$yhat)[1])){  # Loop over MCMC iterations
          yhat_new[m, , ] <- t(sapply(seq_len(dim(bassfunc$yhat)[2]), function(ii) {
            approx(bassfunc$xx.func, bassfunc$yhat[m, ii, ], func.use, rule = 2)$y
          }))
        }
        mod$yhat <- yhat_new

        # Modify the s2 variable?
        sd_full <- sd(bassfunc$y - bassfunc$yhat.mean)
        sd_new  <- sd(mod$y - mod$yhat.mean)
        mod$s2 <- mod$s2 * (sd_new/sd_full)^2
      }
      if(!is.null(mod$yhat.mean)){
        mod$yhat.mean_full <- mod$yhat.mean
        yhat.mean_new <- sapply(seq_len(nrow(mod$yhat.mean)), function(ii) approx(mod$xx.func, mod$yhat.mean[ii,], t_val))
        yhat.mean_new <- unlist(yhat.mean_new[2,])
        mod$yhat.mean <- yhat.mean_new
      }
    }

    # Append the modified model
    results[[cnt]] <- mod
    cnt <- cnt + 1
  }

  # Return a single model if func.use has length 1
  if (length(func.use) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}


# Shifts a vector so nas are at the end
shift_nas <- function(xx){
  not_na_ind <- which(!is.na(xx))
  xx[seq_along(not_na_ind)] <- xx[not_na_ind]
  xx[-seq_along(not_na_ind)] <- NA
  return(xx)
}

# Shifts a matrix so that NAs in the first slice are at the end
shift_nas_group <- function(knots_group) {
  # Identify non-NA rows in the first slice
  not_na_ind <- which(!is.na(knots_group[, 1]))

  # Shift all rows across the third dimension together
  knots_group_tmp <- knots_group[not_na_ind, , drop = FALSE]
  # Add NA rows at the end to maintain original dimensions
  n_na <- nrow(knots_group) - length(not_na_ind)
  knots_group_tmp <- rbind(knots_group_tmp, matrix(NA, nrow = n_na, ncol = dim(knots_group)[2]))

  return(knots_group_tmp)
}
