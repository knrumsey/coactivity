#' FRIEDMAN FUNCTION
#' First input is treated as functional
#' Use p=5, so there is one inert variable
f<-function(x, t){
  10*sin(pi*t*x[1])+20*(x[2]-.5)^2+10*x[3]+5*x[4]
}

#===========================================
#        GENERATE DATA
#===========================================
XX <- lhs::randomLHS(1000, 5)
y1 <- apply(XX, 1, f, t=0.5)
xfunc <- seq(0, 1, length.out=30)
yfunc <- t(apply(XX, 1, f, t=xfunc))

#===========================================
#        CASE 1: Univariate BASS
#===========================================
# This case isn't allowed for augmented version

#===========================================
#      CASE 2: Augmented BASS (fixed t)
#===========================================
mod2 <- bass(XX, yfunc, xx.func=xfunc)

#===========================================
#      CASE 3: PCA BASS (fixed t)
#===========================================
mod3 <- bassPCA(XX, yfunc)

#===========================================
#      ESTIMATE K
#===========================================
C2 <- C_aug_bass(mod2)
C3 <- C_aug_bass(mod3)
