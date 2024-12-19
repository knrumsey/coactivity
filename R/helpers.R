myTimestamp <- function(){
  x <- Sys.time()
  paste("#--", format(x, "%b %d %X"), "--#")
}


scale_range <- function(x, r=NULL){
  if (is.null(r))
    r <- range(x)
  if ((r[2] - r[1]) == 0)
    return(x - r[1])
  return((x - r[1])/(r[2] - r[1]))
}
