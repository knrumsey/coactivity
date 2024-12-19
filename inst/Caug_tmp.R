f <- function(x, t){
  0*t + x[1]^2 + x[1]*x[2]
}
g <- function(x, t){
  u <- 2*t-1
  res <- (f(x, t) + 10*u*x[2]^3) * t * (1 - t)
}
ff <- function(pars) f(pars[1:2], pars[3])
gg <- function(pars) g(pars[1:2], pars[3])


M <- 1e6
Cg <- Cf <- Cfg <-  matrix(0, nrow=3, ncol=3)
for(m in 1:M){
  grad_f <- concordance::fd_grad(ff, runif(3))
  grad_g <- concordance::fd_grad(gg, runif(3))

  Cf <- Cf + tcrossprod(grad_f) / M
  Cg <- Cg + tcrossprod(grad_g) / M
  Cfg <- Cfg + tcrossprod(grad_f, grad_g) / M
  if(m %% 1000 == 0) print(m)
}

# Fit bass models (Old school way)
XX <- lhs::randomLHS(5000, 3)
yf <- apply(XX, 1, ff)
yg <- apply(XX, 1, gg)

modf <- bass(XX, yf)
modg <- bass(XX, yg)


# Functional version
XX <- lhs::randomLHS(1000, 2)
xfunc <- seq(0, 1, length.out=17)
yfunc_f <- t(apply(XX, 1, f, t=xfunc))
yfunc_g <- t(apply(XX, 1, g, t=xfunc))

mod2_f <- bass(XX, yfunc_f, xx.func=xfunc)
mod2_g <- bass(XX, yfunc_g, xx.func=xfunc)


## PCA Version
mod3_f <- bassPCA(XX, yfunc_f)
mod3_g <- bassPCA(XX, yfunc_g)

par(mfrow=c(2,2))
image(Cf, main="Monte Carlo")
image(C_bass(modf), main="coactivity")
image(C_aug_bass(mod2_f), main="functional")
image(C_aug_bass(mod3_f), main="functional")

image(Cg, main="Monte Carlo")
image(C_bass(modg), main="coactivity")
image(C_aug_bass(mod2_g), main="functional")
image(C_aug_bass(mod3_g), main="functional")

image(Cfg, main="Monte Carlo")
image(Cfg_bass(modf, modg), main="coactivity")
image(Cfg_aug_bass(mod2_f, mod2_g))
image(Cfg_aug_bass(mod3_f, mod3_g))


C2 = Cfg_aug_bass(mod2_f, mod2_g)
C3 = Cfg_aug_bass(mod3_f, mod3_g)

eigen(C2)$values / eigen(C3)$values

