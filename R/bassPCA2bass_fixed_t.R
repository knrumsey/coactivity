#' Extract scalar BASS model from bassPCA model for fixed functional variable values
#'
#' Extracts the BASS model(s) corresponding to fixed values of the functional variable.
#' For each specified \( t \)-value, calculates the contribution of the functional variable
#' to the tensor product basis functions and returns a modified BASS model.
#'
#' @param bassPCA A bassBasis model object (from bassPCA() function).
#' @param func.use A numeric vector of fixed values for the functional variable \( t \).
#' @param func.true An optional vector of values for the functional variable in bassPCA. Should have length equal to nrow(bassPCA$dat$basis).
#' @return A list of BASS models for each \( t \)-value (or a single model if length(func.use) == 1).
#' @details
#' Since bassPCA doesn't
#'
#' @examples
#' f<-function(x, t){
#' 10*sin(pi*t*x[1])+20*(x[2]-.5)^2+10*x[3]+5*x[4]
#' }
#' XX <- lhs::randomLHS(500, 5)
#' y1 <- apply(XX, 1, f, t=0.5)
#' xfunc <- seq(0, 1, length.out=20)
#' yfunc <- t(apply(XX, 1, f, t=xfunc))
#'
#' # Fit a bassPCA model
#' mod3_full <- bassPCA(XX, yfunc)
#'
#' # Extract univariate model at t = 0.5
#' mod3 <- bassPCA2bass_fixed_t(mod3_full, 0.5)
#' @export
bassPCA2bass_fixed_t <- function(bassPCA, func.use, func.true=NULL) {
  if (class(bassPCA) != "bassBasis") {
    stop("bassPCA must have class 'bassBasis'.")
  }

  # Extract key quantities
  mod_list <- bassPCA$mod.list
  K <- bassPCA$dat$n.pc
  U <- bassPCA$dat$basis
  Y <- bassPCA$dat$y
  nt <- nrow(U)

  mu <- approx(seq(0,1,length.out=nt),
                 bassPCA$dat$y.m,
                 func.use)$y

  sigma <- approx(seq(0,1,length.out=nt),
                  bassPCA$dat$y.s,
                  func.use)$y

  # Get weights matrix
  wts <- matrix(NA, nrow=length(func.use), ncol=K)
  for(k in 1:K){
    wts[,k] <- get_u(func.use, U, k)
  }

  # Get univariate BASS model
  res <- list()
  cnt <- 1
  for(t_val in func.use){
    # Interpolate to find y, y.m, and y.s
    yy <- apply(Y, 1, function(y) approx(seq(0, 1, length.out=nt),
                                  y,
                                  t_val)$y)

    # Convert linear combination to bass
    res[[cnt]] <- basslc2bass(mod_list,
                              weights=sigma[cnt] * wts[cnt,],
                              offset=mu[cnt],
                              yy=yy)
    cnt <- cnt + 1
  }

  # Return a single model if func.use has length 1
  if (length(func.use) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}
