


colMeans_safe <- function(X) {
  X <- as.matrix(X)
  if(ncol(X)==1){
    return(mean(as.vector(X)))
  }
  return(colMeans(X))
}
bound <- function(x, b){
  pmax(pmin(x,1-b),b)
}
