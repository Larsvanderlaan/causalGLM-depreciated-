
#' Simulate a dataset with known constant CATE
#' @export
#' @param n Sample size
#' @param p Dimension of W
sim.CATE <- function(n=1500, p=2) {
  W <- as.matrix(replicate(p, runif(n, -1,1)))
  colnames(W) <- paste0("W", 1:p)
  betaA <- rnorm(p)
  betaA <- 1.5 * betaA / sum(betaA)
  betaY <- rnorm(p)
  betaY <- 2 * betaY / sum(betaY)
  names(betaA) <- colnames(W)
  names(betaY) <- colnames(betaY)
  g1 <- plogis( W %*% betaA)
  A <- rbinom(n, size = 1, prob = g1)
  
  Q0 <- 0.5*plogis( W %*% betaY)
  
  Q <- 0.45* A  + Q0
  CATE <- 0.45
  
  Y <- rbinom(n, size = 1, prob = Q) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q )
  return(list(descr = "Data simulated from parametric linear model with constant conditional average treatment effect", CATE = CATE, data = data, W = W, A = A, Y= Y, beta_A  = betaA, beta_Y0 = betaY, link = "logistic"))
  
}
#' Simulate a dataset with known constant RR
#' @export
#' @param n Sample size
#' @param p Dimension of W
sim.RR <- function(n=1500, p=2) {
  W <- as.matrix(replicate(p, runif(n, -1,1)))
  colnames(W) <- paste0("W", 1:p)
  betaA <- rnorm(p)
  betaA <- 1.5 * betaA / sum(betaA)
  betaY <- rnorm(p)
  betaY <- 1 * betaY / sum(betaY)
  names(betaA) <- colnames(W)
  names(betaY) <- colnames(betaY)
  
  g1 <- plogis( W %*% betaA)
  A <- rbinom(n, size = 1, prob = g1)
  
  
  Q0 <- 0.1 + 0.3*plogis( W %*% betaY)
  RR <- 2
  Q <- (1-A)*Q0 + A*Q0*RR
  Y <- rbinom(n, size = 1, prob = Q) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q )
  return(list(descr = "Data simulated from parametric linear model with constant relative risk", logRR = log(RR), data = data, W = W, A = A, Y= Y, beta_A  = betaA, beta_Y0 = betaY, link = "logistic"))
}
#' Simulate a dataset with known constant OR
#' @export
#' @param n Sample size
#' @param p Dimension of W
sim.OR <- function(n=1500, p=2) {
  W <- as.matrix(replicate(p, runif(n, -1,1)))
  colnames(W) <- paste0("W", 1:p)
  betaA <- rnorm(p)
  betaA <- 1.5 * betaA / sum(betaA)
  betaY <- rnorm(p)
  betaY <- 2 * betaY / sum(betaY)
  names(betaA) <- colnames(W)
  names(betaY) <- colnames(betaY)
  g1 <- plogis( W %*% betaA)
  A <- rbinom(n, size = 1, prob = g1)
  
  Q0 <- 0.5*plogis( W %*% betaY)
  
  
  Q <- plogis( A - 0.5 +  W %*% betaY)
  OR <- 1
  
  Y <- rbinom(n, size = 1, prob = Q) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q )
  return(list(descr = "Data simulated from parametric linear model with constant odds ratio", logOR = OR, data = data, W = W, A = A, Y= Y, beta_A  = betaA, beta_Y0 = betaY, link = "logistic"))
  
}


 