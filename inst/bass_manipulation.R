# Simulate bass models for each of the three cases
f<-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}

# Functional response (augmented)
n<-500 # number of observations
nfunc<-50 # size of functional variable grid
xfunc<-seq(0,1,length.out=nfunc) # functional grid
x<-matrix(runif(n*9),n,9) # 9 non-functional variables, only first 4 matter
X<-cbind(rep(xfunc,each=n),kronecker(rep(1,nfunc),x)) # to get y
y<-matrix(f(X),nrow=n)

# Fit functional BASS model
mod_t <- bass(x, y, xx.func=xfunc)

# Augment t
mod_aug <- bassfunc2bass(mod_t)
plot(mod_aug)

# Fix t
mod_fix_t <- bassfunc2bass_fixed_t(mod_t, func.use=c(0.5, 0.75))
plot(mod_fix_t[[1]])

# Fit pca model
modPCA_t <- bassPCA(x, y)

# Fix at t
modPCA_fix_t <- bassPCA2bass_fixed_t(modPCA_t, func.use=0.5)



