test_that("simple polynomial example with Gaussian measure", {
  cat('simple polynomial example test with Gaussian measure')

  f <- function(x){
    x[1]^2 + x[1]*x[2] + x[2]^3/9
  }

  C_mc <- function(f, measure, grad=FALSE, nmc=1e4, seed=NULL, ...){
    if(is.null(measure)){
      warning("measure not specified. Setting measure = 5, see help files for details")
      measure <- 5
    }
    if(!is.null(seed)){
      set.seed(seed)
    }
    if(is.numeric(measure)){
      n <- measure[1]
      measure <- function() runif(n)
    }
    if(grad){
      grad_f <- f
    }else{
      grad_f <- function(x, ...) fd_grad(f, x, ...)
    }
    n <- length(measure())

    Cf <- matrix(0, nrow=n, ncol=n)
    for(m in 1:nmc){
      x_m <- measure()
      del_f <- matrix(grad_f(x_m), nrow=n)
      Cf  <- Cf  +  tcrossprod(del_f, del_f)/nmc
    }
    return(Cf)
  }


  # Sim Study Parameters
  N <- 5000
  p <- 3
  sd0 <- 0.05

  # Get true value of C matrix (with monte carlo)
  measure <- function() rnorm(p, 0.5, sd0)
  Cmc <- C_mc(f, measure, nmc=1e5)

  X <- matrix(rnorm(N*p, 0.5, sd0), nrow=N, ncol=p)
  X <- lhs::randomLHS(N, p)
  Yf <- apply(X, 1, f)

  mod_bass <- BASS::bass(X, Yf, verbose=FALSE)
  pr <- list()
  pr[[1]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)
  pr[[2]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)
  pr[[3]]  <- list(dist="normal", trunc=c(-Inf, Inf), mean=0.5, sd=sd0, weights=1)

  Cbass <- C_bass(mod_bass, prior=pr, use_native_scale=TRUE)

  d1 <- sum(abs(Cbass - Cmc))
  expect_that(d1, is_less_than(0.04))
})
