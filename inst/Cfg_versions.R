#' #' FRIEDMAN FUNCTION
#' First input is treated as functional
#' Use p=5, so there is one inert variable
f<-function(x, t){
  10*sin(pi*t*x[1])+20*(x[2]-.5)^2+10*x[3]+5*x[4]
}
g <- function(x, t) {
  5 * sin(2 * pi * t * x[1]) + 15 * (x[3] - 0.3)^2 + 8 * x[4] + 3 * x[5]
}

#===========================================
#        GENERATE DATA
#===========================================
XX <- lhs::randomLHS(500, 5)
y1 <- apply(XX, 1, f, t=0.5)
y2 <- apply(XX, 1, g, t=0.5)

xfunc <- seq(0, 1, length.out=20)
yfunc1 <- t(apply(XX, 1, f, t=xfunc))
yfunc2 <- t(apply(XX, 1, g, t=xfunc))

#===========================================
#        CASE 1: Univariate BASS
#===========================================
mod11 <- bass(XX, y1)
mod12 <- bass(XX, y2)
Cfg1 <- Cfg_bass(mod11, mod12)

#===========================================
#      CASE 2: Augmented BASS (fixed t)
#===========================================
mod2_full1 <- bass(XX, yfunc1, xx.func=xfunc)
mod2_full2 <- bass(XX, yfunc2, xx.func=xfunc)
Cfg2 <- Cfg_bass(mod2_full1, mod2_full2, func.use=0.5)

#===========================================
#      CASE 3: PCA BASS (fixed t)
#===========================================
mod3_full1 <- bassPCA(XX, yfunc1)
mod3_full2 <- bassPCA(XX, yfunc2)
Cfg3 <- Cfg_bass(mod3_full1, mod3_full2, func.use=0.5)

#===========================================
#      CASE 4: COMBINATIONS
#===========================================
Cfg4 <- Cfg_bass(mod11, mod3_full2, func.use=0.5)
Cfg5 <- Cfg_bass(mod11, mod2_full2, func.use=0.5)
Cfg6 <- Cfg_bass(mod2_full1, mod3_full2, func.use=0.5)

#===========================================
#      ESTIMATE Cfg
#===========================================
par(mfrow=c(2,3))
for(i in 1:6){
  image(get(paste0("Cfg",i)))
}
par(mfrow=c(1,1))

# Check that all approaches gives the same eigenvalues
for(i in 1:6){
  if(i == 1){
    plot(eigen(get(paste0("Cfg",i)))$vectors[,1])
  }else{
    points(eigen(get(paste0("Cfg",i)))$vectors[,1], pch=i)
  }
}

# Alternatively
#C2b <- C_bass(mod2_full, func.use = 0.5)
#C3b <- C_bass(mod3_full, func.use = 0.5)

