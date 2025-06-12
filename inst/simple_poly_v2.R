library(BASS)
library(lhs)
#library(coactivity)

ft <- function(x, t) x[1]^2 + x[1]*x[2] + t
gt <- function(x, t) x[1]^2 + x[1]*x[2]*t + t^2

f <- function(x) x[1]^2 + x[1]*x[2] + x[3]
g <- function(x) x[1]^2 + x[1]*x[2]*x[3] + x[3]^2

# Version 1: Treat t as input (manually)
X <- lhs::randomLHS(2000, 3)
yf <- apply(X, 1, f)
yg <- apply(X, 1, g)

modf1 <- bass(X, yf)
modg1 <- bass(X, yg)

Cff <- C_bass(modf1)
Cgg <- C_bass(modg1)
Cfg <- Cfg_bass(modf1, modg1)

# These match perfectly
round(72*Cff)
round(72*Cgg)
round(72*Cfg)

# Version 2: Fit functional bass model
n <- 1000
nfunc <- 100
X <- lhs::randomLHS(n, 2)
t_vec <- seq(0, 1, length.out=nfunc)
XX <- cbind(rep(t_vec, each=n), kronecker(rep(1, nfunc), X))

yft <- t(apply(X, 1, ft, t=t_vec))
ygt <- t(apply(X, 1, gt, t=t_vec))

modf2 <- bass(X, yft, xx.func=t_vec)
modg2 <- bass(X, ygt, xx.func=t_vec)

round(C_aug_bass(modf2)*72)
round(C_aug_bass(modg2)*72)
round(Cfg_aug_bass(modf2, modg2)*72)

# Version 3: Fit PCA bass model
modf3 <- bassPCA(X, yft, center=FALSE, n.pc=2)
modg3 <- bassPCA(X, ygt, center=FALSE, n.pc=3)

round(72*100*C_aug_bass(modf3))
round(72*100*C_aug_bass(modg3))
round(72*100*Cfg_aug_bass(modf3, modg3))



