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


xfunc <- seq(0, 1, length.out=10)
yfunc1 <- t(apply(XX, 1, f, t=xfunc))
yfunc2 <- t(apply(XX, 1, g, t=xfunc))

#===========================================
#        CASE 1: Univariate BASS
#===========================================
mod1_list <- mod2_list <- list()
cnt <- 1
for(t in xfunc){
  yt1 <- apply(XX, 1, f, t=t)
  yt2 <- apply(XX, 1, g, t=t)
  mod1_list[[cnt]] <- bass(XX, yt1)
  mod2_list[[cnt]] <- bass(XX, yt2)
  cnt <- cnt + 1
}

#===========================================
#      CASE 2: Augmented BASS (fixed t)
#===========================================
mod2_full1 <- bass(XX, yfunc1, xx.func=xfunc)
mod2_full2 <- bass(XX, yfunc2, xx.func=xfunc)

#===========================================
#      CASE 3: PCA BASS (fixed t)
#===========================================
mod3_full1 <- bassPCA(XX, yfunc1)
mod3_full2 <- bassPCA(XX, yfunc2)

#===========================================
#      ESTIMATE Kfg
#===========================================
Kfg_bass(mod1_list, mod2_list, func.use=xfunc)
Kfg_bass(mod2_full1, mod2_full1, func.use=xfunc)
Kfg_bass(mod3_full1, mod3_full2, func.use=xfunc)
# Combinations
Kfg_bass(mod1_list,  mod3_full2, func.use=xfunc)
Kfg_bass(mod2_full1, mod3_full2, func.use=xfunc)
