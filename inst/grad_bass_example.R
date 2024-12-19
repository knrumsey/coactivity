# dms_additive
f <- function(x){
  1.3356 * (1.5 * (1 - x[1]) + exp(2 * x[1] - 1) * sin(3 * pi *
           (x[1] - 0.6)^2) + exp(3 * (x[2] - 0.5)) * sin(4 * pi * (x[2] - 0.9)^2))
}
# Gradient of dms_additive
gradient_f <- function(x) {
  x1 <- x[1]
  x2 <- x[2]

  # Partial derivative with respect to x1
  term1 <- -1.5
  term2 <- 2 * exp(2 * x1 - 1) * sin(3 * pi * (x1 - 0.6)^2)
  term3 <- exp(2 * x1 - 1) * cos(3 * pi * (x1 - 0.6)^2) * 6 * pi * (x1 - 0.6)
  grad_x1 <- 1.3356 * (term1 + term2 + term3)

  # Partial derivative with respect to x2
  term4 <- 3 * exp(3 * (x2 - 0.5)) * sin(4 * pi * (x2 - 0.9)^2)
  term5 <- exp(3 * (x2 - 0.5)) * cos(4 * pi * (x2 - 0.9)^2) * 8 * pi * (x2 - 0.9)
  grad_x2 <- 1.3356 * (term4 + term5)

  # Return gradient as a vector
  c(grad_x1, grad_x2)
}


X <- lhs::randomLHS(500, 2)
y <- apply(X, 1, f)
fhat <- BASS::bass(X, y)

grads <- gradient_bass(fhat, X, mcmc.use = seq(1, 1000, by=50), verbose=TRUE)
grads_bass <- t(apply(grads, c(2,3), mean))
grads_true <- apply(X, 1, gradient_f)
plot(unlist(grads_bass), unlist(grads_true))
abline(0, 1, lwd=2, col='orange')


# Super simple example
pos <- BASS:::pos
g <- function(x){
  5*pos(x[1]-0.5) - 2*pos(x[2] - 0.7)
}
gradient_g <- function(x) {
  grad_x1 <- ifelse(x[1] > 0.5, 5, 0)
  grad_x2 <- ifelse(x[2] > 0.7, -2, 0)
  c(grad_x1, grad_x2)
}

X <- lhs::randomLHS(500, 2)
y <- apply(X, 1, g)
ghat <- BASS::bass(X, y, maxBasis = 2, temp.ladder = 1.02^(0:10))

grads <- gradient_bass(ghat, X, mcmc.use = seq(1, 1000, by=50), verbose=TRUE)
grads_bass <- t(apply(grads, c(2,3), mean))
grads_true <- apply(X, 1, gradient_g)
plot(as.numeric(grads_bass), as.numeric(grads_true))
abline(0, 1, lwd=2, col='orange')


