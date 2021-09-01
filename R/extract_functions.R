#' @export
summary.causalGLM <- function(object) {
  (object$coefs[,c(1,6,4,5,7)])
}

#' @export
coef.causalGLM<- function(object) {
  (object$coefs)
}


#' @export
predict.causalGLM <- function(object, Wnew) {
  (object$pred_function(Wnew))
}
