library(BASS)
library(lhs)
#library(coactivity)

ft <- function(x, t) x[1] * t
gt <- function(x, t) 2 * x[1] * t^2

n <- 500
nfunc <- 100
X <- lhs::randomLHS(n, 1)
t_vec <- seq(0, 1, length.out=nfunc)
XX <- cbind(rep(t_vec, each=n), kronecker(rep(1, nfunc), X))

yft <- t(apply(X, 1, ft, t=t_vec)) + rnorm(n*nfunc, 0, 1e-4)
ygt <- t(apply(X, 1, gt, t=t_vec)) + rnorm(n*nfunc, 0, 1e-4)

# Version 2
modf2 <- bass(X, yft, xx.func=t_vec)
modg2 <- bass(X, ygt, xx.func=t_vec)

round(C_aug_bass(modf2)*6)
round(C_aug_bass(modg2)*6)
round(Cfg_aug_bass(modf2, modg2)*6)

# Version 3: Fit PCA bass model
modf3 <- bassPCA(X, yft, center=FALSE)
modg3 <- bassPCA(X, ygt, center=FALSE)

round(6*100*C_aug_bass(modf3))
round(6*100*C_aug_bass(modg3))
round(6*100*Cfg_aug_bass(modf3, modg3))

