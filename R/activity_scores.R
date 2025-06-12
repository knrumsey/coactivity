#' Activity Scores
#'
#' This function computes the activity scores for main effects of the variables
#'
#' @param obj Either a BASS model, a matrix, or a list of matrices (with one element per MCMC iteration)
#' @param k The number of columns of W to consider
#' @param mcmc.use a vector of indices telling which MCMC draws to use
#' @param norm Logical, should activity scores be normalized to have length one?
#' @param plot_it Logical, should a plot be made?
#' @param ... additional arguments passed to the corresponding `C_bass` function.
#' @return the activity scores
#' @export
activity_scores <- function(obj, k = 1,  mcmc.use = NULL, norm = FALSE, plot_it = FALSE, ...) {
  # Handle case: BASS model object
  if("bass" %in% class(obj) || "bassBasis" %in% class(obj)){
    # Get list of C matrices
    C_list <- C_bass(obj, mcmc.use = mcmc.use, ...)

    # If C_list is a single matrix (not list), coerce into list
    if(is.matrix(C_list)) {
      C_list <- list(C_list)
    }
  }else if(is.matrix(obj)){
    C_list <- list(obj)  # Single matrix as input
  }else if (is.list(obj) && all(sapply(obj, is.matrix))){
    if(is.null(mcmc.use)){
      C_list <- obj
    }else{
      C_list <- obj[mcmc.use]
    }
  }else{
    stop("Input 'obj' must be a BASS model, a matrix, or a list of matrices.")
  }

  # Compute activity scores for each matrix
  score_list <- lapply(C_list, function(C){
    eig <- eigen(C, symmetric=TRUE)
    W <- eig$vectors
    lam <- eig$values

    # Calculate activity scores
    Lam_k <- diag(pmax(lam[1:k], 0), k)
    W_k <- W[, 1:k, drop=FALSE]
    res <- rowSums(W_k^2 %*% Lam_k)

    # normalize if requested
    if(norm){
      res <- res / sqrt(sum(res^2))
    }
    return(res)
  })

  # Reduce if only one element
  if(length(score_list) == 1) {
    score_list <- score_list[[1]]
  }else{
    score_list <- do.call(rbind, score_list)
  }

  # Plot if requested (only for single result or averaged)
  if(is.matrix(score_list)){
    #browser()
    p <- ncol(score_list)
    n <- nrow(score_list)

    plot(NULL,
         xlim=c(1-0.2, p+0.2), ylim=range(score_list),
         xaxt="n",
         ylab=paste0("Activity Score (", k, ")"), xlab="Input")
    axis(1, 1:p, 1:p)
    for(i in 1:p){
      points(rnorm(n, i, 0.025), score_list[,i], col="grey")
      points(i, mean(score_list[,i]), pch=17)
    }
  }else{
    # Single score vector
    plot(score_list,
         xlab = "Inputs", ylab = paste0("Activity Score (", k, ")"),
         pch = 16, cex = 2, col = "black",
         ylim = c(0, max(score_list)))
  }

  return(score_list)
}



#' Co-Activity Scores
#'
#' This function computes the co-activity scores for main effects of the variables
#'
#' @param obj1 Either a BASS model, a matrix, or a list of matrices (with one element per MCMC iteration).
#' @param obj2 Either a BASS model or NULL. If NULL, then \code{obj1} should be a matrix (or list of matrices) corresponding to the Vfg matrix.
#' @param k The number of columns of W to consider
#' @param signed Use signed or unsigned version of coactivity scores?
#' @param mcmc.use a vector of indices telling which MCMC draws to use
#' @param norm Logical, should activity scores be normalized to have length one?
#' @param plot_it Logical, should a plot be made?
#' @param ... additional arguments passed to the corresponding `Cfg_bass` function.
#' @return the co-activity scores
#' @export
coactivity_scores <- function(obj1, obj2 = NULL, k = 1, signed = TRUE, mcmc.use = NULL, norm = FALSE, plot_it = TRUE, ...) {
  # Step 1: Get the co-Constantine matrix/matrices (Vfg)
  if (!is.null(obj2)) {
    # Case: both obj1 and obj2 are BASS models
    Vfg_list <- Cfg_bass(obj1, obj2, mcmc.use = mcmc.use, ...)

    if(is.matrix(Vfg_list)){
      Vfg_list <- list(Vfg_list)
    }

    for(i in seq_along(Vfg_list)){
      VV <- Vfg_list[[i]]
      Vfg_list[[i]] <- VV + t(VV)
    }

  }else{
    # Case: Vfg is supplied directly as matrix or list of matrices
    if(is.matrix(obj1)){
      Vfg_list <- list(obj1)
    }else if(is.list(obj1) && all(sapply(obj1, is.matrix))){
      Vfg_list <- obj1
    }else{
      stop("If obj2 is NULL, then obj1 must be a matrix or list of matrices.")
    }
  }

  # Step 2: Compute co-activity scores for each matrix
  score_list <- lapply(Vfg_list, function(V) {
    eig <- eigen(V, symmetric = TRUE)
    W <- eig$vectors
    lam <- eig$values

    if(!signed){
      lam <- abs(lam)
    }

    ord <- rev(order(abs(lam)))
    lam_ord <- lam[ord]
    W_ord <- W[, ord, drop = FALSE]

    Lam_k <- diag(lam_ord[1:k], k)
    W_k <- W_ord[, 1:k, drop = FALSE]
    res <- rowSums(W_k^2 %*% Lam_k)

    if(norm){
      res <- res / sqrt(sum(res^2))
    }

    return(res)
  })

  # Step 3: Reduce if only one matrix
  if (length(score_list) == 1) {
    score_list <- score_list[[1]]
  } else {
    score_list <- do.call(rbind, score_list)
  }

  # Step 4: Plot if requested
  if (plot_it) {
    if (is.matrix(score_list)) {
      p <- ncol(score_list)
      n <- nrow(score_list)

      plot(NULL,
           xlim = c(1 - 0.2, p + 0.2), ylim = range(score_list),
           xaxt = "n",
           ylab = paste0("Coactivity Score (", k, ")"), xlab = "Input")
      axis(1, 1:p, 1:p)
      for (i in 1:p) {
        points(rnorm(n, i, 0.025), score_list[, i], col = "grey")
        points(i, mean(score_list[, i]), pch = 17)
      }
    } else {
      plot(score_list,
           xlab = "Inputs", ylab = paste0("Coactivity Score (", k, ")"),
           pch = 17, col = "black",
           ylim = c(0, max(score_list)))
    }
  }

  return(score_list)
}


#' Concordance
#'
#' This function computes the concordance between two models
#'
#' @param obj1 Either a BASS model or a matrix (or a list of matrices) with one element per MCMC iteration, corresponding to Cf.
#' @param obj2 Either a BASS model or a matrix (or a list of matrices) with one element per MCMC iteration, corresponding to Cg.
#' @param Cfg A matrix or list of matrices corresponding to Cfg. Required if obj1/obj2 are not both BASS models.
#' @param mcmc.use a vector of indices telling which MCMC draws to use
#' @param ... additional arguments passed to the corresponding `C_bass` or `Cfg_bass` function.
#' @return A scalar or vector of concordance scores
#' @export
concordance <- function(obj1, obj2, Cfg = NULL, mcmc.use = NULL, ...) {
  # Step 1: If both inputs are BASS models
  if (("bass" %in% class(obj1) || "bassBasis" %in% class(obj1)) &&
      ("bass" %in% class(obj2) || "bassBasis" %in% class(obj2))) {

    Cf_list <- C_bass(obj1, mcmc.use = mcmc.use, ...)
    Cg_list <- C_bass(obj2, mcmc.use = mcmc.use, ...)
    Cfg_list <- Cfg_bass(obj1, obj2, mcmc.use = mcmc.use, ...)

    # Ensure each is a list
    if (is.matrix(Cf_list)) Cf_list <- list(Cf_list)
    if (is.matrix(Cg_list)) Cg_list <- list(Cg_list)
    if (is.matrix(Cfg_list)) Cfg_list <- list(Cfg_list)

  } else {
    # Otherwise, we expect all to be provided
    if (is.null(Cfg)) stop("If obj1 and obj2 are not both BASS models, Cfg must be supplied.")

    # Coerce each input to a list of matrices
    to_list <- function(x) if (is.matrix(x)) list(x) else if (is.list(x) && all(sapply(x, is.matrix))) x else stop("Each object must be a matrix or list of matrices.")

    Cf_list <- to_list(obj1)
    Cg_list <- to_list(obj2)
    Cfg_list <- to_list(Cfg)
  }

  # Step 2: Check lengths
  lens <- c(length(Cf_list), length(Cg_list), length(Cfg_list))
  min_len <- min(lens)

  if (length(unique(lens)) != 1) {
    warning("Cf, Cg, and Cfg have different lengths; truncating to minimum length (", min_len, ").")
    Cf_list <- Cf_list[1:min_len]
    Cg_list <- Cg_list[1:min_len]
    Cfg_list <- Cfg_list[1:min_len]
  }

  # Step 3: Compute concordance per iteration
  concordance_vec <- vapply(seq_len(min_len), function(i) {
    Cf <- Cf_list[[i]]
    Cg <- Cg_list[[i]]
    Cfg <- Cfg_list[[i]]

    tr_Cfg <- sum(diag(Cfg))
    tr_Cf <- sum(diag(Cf))
    tr_Cg <- sum(diag(Cg))

    if (tr_Cf <= 0 || tr_Cg <= 0) return(NA_real_)  # Prevent division by zero
    tr_Cfg / sqrt(tr_Cf * tr_Cg)
  }, numeric(1))

  # Step 4: Return scalar or vector
  if (length(concordance_vec) == 1) {
    return(concordance_vec[[1]])
  } else {
    return(concordance_vec)
  }
}

#' Gradient-Based Shapley Values
#'
#' This function computes the Shapley Values (gradient-based) as described by Duan and Okten (2023).
#'
#' @param obj Either a BASS model, a matrix, or a list of matrices (with one element per MCMC iteration)
#' @param mcmc.use a vector of indices telling which MCMC draws to use
#' @param norm Logical, should Shapley values be normalized to have length one?
#' @param plot_it Logical, should a plot be made?
#' @param ... additional arguments passed to the corresponding `C_bass` function.
#' @references Duan, Hui, and Giray Okten. "Derivative-based Shapley value for global sensitivity analysis and machine learning explainability." International Journal for Uncertainty Quantification 15.1 (2025).
#' @return the Shapley values
#' @export
shapley_grad <- function(obj, mcmc.use = NULL, norm = FALSE, plot_it = FALSE, ...) {
  # Coerce to list of C matrices
  if ("bass" %in% class(obj) || "bassBasis" %in% class(obj)) {
    C_list <- C_bass(obj, mcmc.use = mcmc.use, ...)
    if (is.matrix(C_list)) C_list <- list(C_list)
  } else if (is.matrix(obj)) {
    C_list <- list(obj)
  } else if (is.list(obj) && all(sapply(obj, is.matrix))) {
    C_list <- if (is.null(mcmc.use)) obj else obj[mcmc.use]
  } else {
    stop("Input 'obj' must be a BASS model, a matrix, or a list of matrices.")
  }

  # Compute Shapley values per matrix
  score_list <- lapply(C_list, function(C) {
    absC <- abs(C)
    phi <- rowSums(absC)
    diagC <- diag(C)
    phi <- 0.5 * (phi + diagC)
    if (norm) {
      phi <- phi / sqrt(sum(phi^2))
    }
    return(phi)
  })

  # Reduce
  if (length(score_list) == 1) {
    score_list <- score_list[[1]]
  } else {
    score_list <- do.call(rbind, score_list)
  }

  # Plot if requested
  if (plot_it) {
    if (is.matrix(score_list)) {
      p <- ncol(score_list)
      n <- nrow(score_list)

      plot(NULL,
           xlim = c(1 - 0.2, p + 0.2), ylim = range(score_list),
           xaxt = "n",
           ylab = "Shapley Value", xlab = "Input")
      axis(1, 1:p, 1:p)
      for (i in 1:p) {
        points(rnorm(n, i, 0.025), score_list[, i], col = "grey")
        points(i, mean(score_list[, i]), pch = 17)
      }
    } else {
      plot(score_list,
           xlab = "Inputs", ylab = "Shapley Value",
           pch = 16, col = "black",
           ylim = c(0, max(score_list)))
    }
  }

  return(score_list)
}


#' Gradient-Based Co-Shapley Values
#'
#' This function computes co-Shapley values from the cross-partial covariance matrix Cfg.
#'
#' @param obj1 A BASS model, or a matrix or list of matrices representing Cfg.
#' @param obj2 Either a BASS model (if computing Cfg from two models), or NULL (if Cfg is passed in directly).
#' @param mcmc.use Vector of MCMC draws to use (for BASS models).
#' @param use_Vfg Logical (default TRUE). Should the symmetrized matrix Vfg = (Cfg + Cgf)/2 be used?
#' @param norm Logical. Should values be normalized to unit length?
#' @param plot_it Logical. Should a plot be produced?
#' @param ... Additional arguments passed to \code{Cfg_bass()} if computing from BASS models.
#' @details The analogous quantity for Cfg rather than Cf.
#' @return A vector or matrix of co-Shapley scores.
#' @export
coshapley_grad <- function(obj1, obj2 = NULL, mcmc.use = NULL, use_Vfg=TRUE, norm = FALSE, plot_it = FALSE, ...) {
  # Step 1: Construct the list of Cfg matrices
  if (!is.null(obj2)) {
    # obj1 and obj2 are BASS models — compute Cfg
    Cfg_list <- Cfg_bass(obj1, obj2, mcmc.use = mcmc.use, ...)
    if (is.matrix(Cfg_list)) Cfg_list <- list(Cfg_list)

    # Symmetrize Cfg
    if(use_Vfg){
      for (i in seq_along(Cfg_list)) {
        Cfg_list[[i]] <- (Cfg_list[[i]] + t(Cfg_list[[i]]))/2
      }
    }
  } else {
    # Assume obj1 is a matrix or list of matrices representing Cfg
    if (is.matrix(obj1)) {
      Cfg_list <- list(obj1)
    } else if (is.list(obj1) && all(sapply(obj1, is.matrix))) {
      Cfg_list <- obj1
      if (!is.null(mcmc.use)) {
        Cfg_list <- Cfg_list[mcmc.use]
      }
    } else {
      stop("If obj2 is NULL, obj1 must be a matrix or list of matrices representing Cfg.")
    }
  }

  # Step 2: Compute co-Shapley values
  score_list <- lapply(Cfg_list, function(Cfg) {
    absCfg <- abs(Cfg)
    phi <- rowSums(absCfg)
    diagCfg <- diag(Cfg)
    phi <- 0.5 * (phi + diagCfg)
    if (norm) {
      phi <- phi / sqrt(sum(phi^2))
    }
    return(phi)
  })

  # Step 3: Reduce if single matrix
  if (length(score_list) == 1) {
    score_list <- score_list[[1]]
  } else {
    score_list <- do.call(rbind, score_list)
  }

  # Step 4: Plot if requested
  if (plot_it) {
    if (is.matrix(score_list)) {
      p <- ncol(score_list)
      n <- nrow(score_list)

      plot(NULL,
           xlim = c(1 - 0.2, p + 0.2), ylim = range(score_list),
           xaxt = "n",
           ylab = "Co-Shapley Value", xlab = "Input")
      axis(1, 1:p, 1:p)
      for (i in 1:p) {
        points(rnorm(n, i, 0.025), score_list[, i], col = "grey")
        points(i, mean(score_list[, i]), pch = 17)
      }
    } else {
      plot(score_list,
           xlab = "Inputs", ylab = "Co-Shapley Value",
           pch = 16, col = "black",
           ylim = c(0, max(score_list)))
    }
  }

  return(score_list)
}
