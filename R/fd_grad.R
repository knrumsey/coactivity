#' Forward auto-differentiation function
#'
#' Function for approximating the gradient of a function
#'
#' @param f the function to find the gradient of
#' @param x the input values
#' @param h the tolerance
#' @param ... additional inputs to be passed to f
#' @return The approximate gradient of f at x
#' @export
fd_grad <- function(f, x, h=1e-12, ...){
  gradient <- rep(NA, length(x))
  for(i in seq_along(x)){
    x[i] <- x[i] + h*1i
    gradient[i] <- Im(f(x, ...))/h
    x[i] <- x[i] - h*1i
  }
  return(gradient)
}
