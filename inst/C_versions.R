#' FRIEDMAN FUNCTION
#' First input is treated as functional
#' Use p=5, so there is one inert variable
f<-function(x, t){
  10*sin(pi*t*x[1])+20*(x[2]-.5)^2+10*x[3]+5*x[4]
}

#===========================================
#        GENERATE DATA
#===========================================
XX <- lhs::randomLHS(500, 5)
y1 <- apply(XX, 1, f, t=0.5)
xfunc <- seq(0, 1, length.out=20)
yfunc <- t(apply(XX, 1, f, t=xfunc))

Xtest <- lhs::randomLHS(100, 5)
ytest <- apply(Xtest, 1, f, t=0.5)

#===========================================
#        CASE 1: Univariate BASS
#===========================================
mod1 <- bass(XX, y1)

#===========================================
#      CASE 2: Augmented BASS (fixed t)
#===========================================
mod2_full <- bass(XX, yfunc, xx.func=xfunc)
mod2 <- bassfunc2bass_fixed_t(mod2_full, 0.5)

#===========================================
#      CASE 3: PCA BASS (fixed t)
#===========================================
mod3_full <- bassPCA(XX, yfunc)
mod3 <- bassPCA2bass_fixed_t(mod3_full, 0.5)

#===========================================
#      VALIDATE MODELS
#===========================================
yhat1 <- colMeans(predict(mod1, Xtest))
yhat2 <- colMeans(predict(mod2, Xtest))
yhat3 <- colMeans(predict(mod3, Xtest, mcmc.use=1:1000))

par(mfrow=c(1,3))
for(i in 1:3){
  plot(ytest, get(paste0("yhat",i)))
  abline(0,1,col='orange')
}

#===========================================
#      ESTIMATE C
#===========================================
C1 <- C_bass(mod1)
C2 <- C_bass(mod2)
C3 <- C_bass(mod3)

image(C1)
image(C2)
image(C3)
par(mfrow=c(1,1))

# Alternatively
#C2b <- C_bass(mod2_full, func.use = 0.5)
#C3b <- C_bass(mod3_full, func.use = 0.5)

