test_that("simple polynomial example with multivariate Gaussian measure", {
  cat('simple polynomial example test with multivariate Gaussian measure')

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
  N <- 1000
  p <- 3

  # Prior
  mu <- rep(0.5, 3)
  Sigma <- matrix(c(1, 0.5, -0.2,
                    0.5, 1, 0.1,
                    -0.2, 0.1, 1),
                  ncol=3, byrow=TRUE)/3

  # MC Estimate
  measure <- function() as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))
  Cmc4 <- C_mc(f, measure, nmc=5e4)

  # BASS Estimate
  pr <- list()
  for(i in 1:3) pr[[i]] <- list(dist="normal", trunc=c(-Inf, Inf), mean=0, sd=1, weights=1)
  X <- lhs::randomLHS(N, p)*4.5 - 2
  y <- apply(X, 1, f)

  Esig <- eigen(Sigma)
  A <- Esig$vectors%*%diag(1/sqrt(Esig$values))%*%t(Esig$vectors)
  X0 <- X - rep(mu, each=N)
  Z <- X0%*%A

  mod4 <- BASS::bass(Z, y, nmcmc=25000, nburn=20000, thin=5)
  Cba4z <- C_bass(mod4, pr)
  Cba4 <- t(A)%*%Cba4z%*%A

  d1 <- max(abs(Cba4-Cmc4))/max(Cba4)
  expect_that(d1, is_less_than(0.15))
})
