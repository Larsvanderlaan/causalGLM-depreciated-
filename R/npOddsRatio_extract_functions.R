#' @export
summary.npOddsRatio <- function(object) {
  (object$coefs[,c(1,6,4,5,7)])
}


coef.npOddsRatio <- function(object) {
  (object$coefs)
}


#' @export
predict.npOddsRatio <- function(object, Wnew) {
  (object$pred_function(Wnew))
}
