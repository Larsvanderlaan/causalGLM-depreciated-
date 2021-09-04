
 #    out <-  sim.causalGLM(cross_fit = F, formula = ~1 + W1 + W2 + W3 + W4, n=50, p = 4, learning_method = "glmnet", nsims = 500, estimand= "CATE", compare_with_competitor = T)
#' @export
sim.causalGLM <- function(formula = ~ 1 + W1 + W2 + W3 + W4 ,  estimand = c("CATE", "OR", "RR"), n = 75, p = 4, prop_active = 1, nsims = 100, learning_method = "glmnet", formula_true = formula,  formula_A = ~., formula_Y0W = ~., silent = FALSE, compare_with_competitor = TRUE, print_CI_length = F  , ... ){
 
  estimand <- match.arg(estimand)
  if(estimand == "CATE") {
    sim <- sim.CATE
  } else if (estimand == "OR") {
    sim <- sim.OR
  } else if (estimand == "RR") {
    sim <- sim.RR
  }
  data_list <- sim(n,p,prop_active,  formula_estimand = formula_true, formula_A = formula_A, formula_Y0W = formula_Y0W)
  beta <- data_list[[2]]
  print(beta)
  beta_A <- data_list$beta_A
  beta_Y <- data_list$beta_Y
  passes <- c()
  passes1 <- c()
  #type1 <- c()
  print(paste0("Estimand is ", estimand))
  print(paste0("n=", n, "; p=",p))
  print(paste0("Total sims is: ", nsims))
   
  for(i in 1:nsims) {
    print(paste0("Running simulation iteration: ", i))
    data_list <- sim(n,p,prop_active,  formula_estimand = formula_true, formula_A = formula_A, formula_Y0W = formula_Y0W, beta = beta, beta_A = beta_A, beta_Y=beta_Y)
    
      fit_glm <- causalGLM(formula, W = data_list$W, A = data_list$A, Y = data_list$Y, estimand = estimand, learning_method = learning_method, constant_variance_CATE = TRUE , return_competitor = compare_with_competitor, ...)
      
      
    
     out <- coef(fit_glm)
     
    if(is.vector(out)) {
      matrix(out, nrow=1)
    }
    if(i%%25==0){
     # print(out)
    }
    ci <- out[,c(4,5), drop =F]
    if(print_CI_length) {
      print("causlglm")
      print(out[,1])
      print(ci)
      print(apply(ci,1,diff))
    }
    pval <-  out[,c(7)]
    passes <- cbind(passes, ci[,1] <= beta & ci[,2] >= beta )
     
     
    print("Coverage probability of 95% confidence intervals of causalGLM so far: ")
    print(rowMeans(passes))
    
    out <- fit_glm$coefs1
    if(!is.null(out)) {
       
      ci <- out[,c(4,5), drop =F]
      if(print_CI_length) {
        print("est")
        print(out[,1])
        print(ci)
        print(apply(ci,1,diff))
      }
      passes1 <- cbind(passes1, ci[,1] <= beta & ci[,2] >= beta )
        
        
      print("Coverage probability of 95% confidence intervals of the estimating equation competitor so far: ")
      print(rowMeans(passes1))
     
      
    }
     
    # print("Proportion of p-values less than 0.05 so far: ")
    # print(mean(type1))
    report <- function() {
      print("Coverage probability of 95% confidence intervals of causalGLM so far: ")
      print(rowMeans(passes))
      print("Coverage probability of 95% confidence intervals of the estimating equation competitor so far: ")
      print(rowMeans(passes1))
    }
  }
  return(list(report = report, CI_coverage = rowMeans(passes), CI_coverage_competitor= rowMeans(passes1), alpha = 0.05, data_list = data_list))
   
  
}

#     out <-   sim.causalGLMwithLasso(compare_with_competitor = T, nsims = 1000, cross_fit = TRUE)
bigformula <- ~ 1 + W12 + + W2 + W3 + W112 + W123 + W254 + W345 + W234 + W421  + W163 + W164 + W262 + W263 + W419 +  W22 + W3 + W4 + W50 + W86 + W25 + W118 + W299 + W10 + W110 + W312 + W213 + W414 + W290+ W16 + W150 + W333 + W222 + W456
smallformula <- ~ 1 + W12 + W50  + W299 + W414 
#' @export
sim.causalGLMwithLasso <- function(formula = smallformula,  estimand = c("CATE", "OR", "RR"), n = 100, p = 500, prop_active = 1, nsims = 100, learning_method = "glmnet", formula_true = formula,  formula_A = bigformula, formula_Y0W = bigformula, silent = FALSE, compare_with_competitor = F,    ... ){
  
  estimand <- match.arg(estimand)
  if(estimand == "CATE") {
    sim <- sim.CATE
  } else if (estimand == "OR") {
    sim <- sim.OR
  } else if (estimand == "RR") {
    sim <- sim.RR
  }
  data_list <- sim(n,p,prop_active,  formula_estimand = formula_true, formula_A = formula_A, formula_Y0W = formula_Y0W)
  beta <- data_list[[2]]
  print(beta)
  beta_A <- data_list$beta_A
  beta_Y <- data_list$beta_Y
  passes <- c()
  passes1 <- c()
  type1 <- c()
  print(paste0("Estimand is ", estimand))
  print(paste0("n=", n, "; p=",p))
  print(paste0("Total sims is: ", nsims))
  print(dim(data_list$W))
  for(i in 1:nsims) {
    print(paste0("Running simulation iteration: ", i))
    data_list <- sim(n,p,prop_active,  formula_estimand = formula_true, formula_A = formula_A, formula_Y0W = formula_Y0W, beta = beta, beta_A = beta_A, beta_Y=beta_Y)
    
    fit_glm <- causalGLMwithLASSO(formula, W = data_list$W, A = data_list$A, Y = data_list$Y, estimand = estimand,   constant_variance_CATE = TRUE , return_competitor = compare_with_competitor, ...)
    
    
    
    out <- coef(fit_glm)
    if(is.vector(out)) {
      matrix(out, nrow=1)
    }
    if(i%%25==0){
      print(out)
    }
    ci <- out[,c(4,5), drop =F]
    pval <-  out[,c(7)]
    passes <- cbind(passes, ci[,1] <= beta & ci[,2] >= beta )
    out <- fit_glm$coefs1
    
    print("Coverage probability of 95% confidence intervals of causalGLM so far: ")
    print(rowMeans(passes))
    if(!is.null(out)) {
      ci <- out[,c(4,5), drop =F]
      passes1 <- cbind(passes1, ci[,1] <= beta & ci[,2] >= beta )
      print("Coverage probability of 95% confidence intervals of the estimating equation  competitor so far: ")
      print(rowMeans(passes1))
    }
    
    # print("Proportion of p-values less than 0.05 so far: ")
    # print(mean(type1))
  }
  report <- function() {
    print("Coverage probability of 95% confidence intervals of causalGLM so far: ")
    print(rowMeans(passes))
    print("Coverage probability of 95% confidence intervals of the estimating equation (DML) competitor so far: ")
    print(rowMeans(passes1))
  }
  return(list(report = report, CI_coverage = rowMeans(passes),   alpha = 0.05, data_list = data_list))
  
  
}
  
  
#' Simulate a dataset with known constant CATE
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.CATE <- function(n=1500, p=2, prop_active = 1,  sigma = NULL, formula_estimand = ~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -1,1)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A,as.data.frame(W))
  XY <- model.matrix(formula_Y0W,as.data.frame(W))
  V <- model.matrix(formula_estimand,as.data.frame(W))
  pA<- ncol(XA)
  pY <- ncol(XY)
   
  if(is.null(beta_A)) {
  activeA <- rbinom(pA, size = 1, prob =prop_active)
  beta_A <- runif(pA, min=-1,max=1)
  beta_A <- 1.5 * beta_A * activeA / sum(abs(beta_A*activeA))
  names(beta_A) <- colnames(XA)
  }
  if(is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob =prop_active)
    beta_Y <- runif(pY, min=-1,max=1)
    beta_Y <- 2 * beta_Y*activeY / sum(abs(beta_Y*activeY))
    
    names(beta_Y) <- colnames(XY)
  }
  
  g1 <- plogis( XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)
  Q0 <-  XY %*% beta_Y
  
  if(!is.null(beta)) {
    beta_CATE <- beta
  } else {
    beta_CATE <- runif(ncol(V), min=-1,max=1)
    beta_CATE <-  2*beta_CATE / sum(abs(beta_CATE))
  }
   
  CATE <- V %*% beta_CATE
  
  Q <-CATE* A  + Q0
  if(is.null(sigma)) {
    sigma <- sd(Q)/4
  }
  Y <- rnorm(n, mean = Q, sd = sigma) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q , sd = sigma,  CATE = CATE)
  return(list(descr = "Data simulated from parametric linear model with constant conditional average treatment effect", beta_CATE = beta_CATE, data = data, W = W, A = A, Y= Y, beta_A  = beta_A, beta_Y0 = beta_Y, link = "logistic"))
  
}
#' Simulate a dataset with known constant RR
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.RR <- function(n=1500, p=2, prop_active = 1, formula_estimand =~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -0.5,0.5)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A,as.data.frame(W))
  XY <- model.matrix(formula_Y0W,as.data.frame(W))
  V <- model.matrix(formula_estimand,as.data.frame(W))
  pA<- ncol(XA)
  pY <- ncol(XY)
  if(is.null(beta_A)) {
    activeA <- rbinom(pA, size = 1, prob =prop_active)
    beta_A <- runif(pA, min=-1,max=1)
    beta_A <- 3 * beta_A * activeA / sum(abs(beta_A*activeA))
    names(beta_A) <- colnames(XA)
  }
  if(is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob =prop_active)
    beta_Y <- runif(pY, min=0,max=1)
    beta_Y <- 0.2 + 3 * beta_Y*activeY / sum(abs(beta_Y*activeY))
    names(beta_Y) <- colnames(XY)
  }
  g1 <- plogis( XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)
  Q0 <- exp( XY %*% beta_Y)
  if(!is.null(beta)) {
    betalogRR <- beta
  } else {
    betalogRR <- runif(ncol(V), min=0,max=1)
    betalogRR <- 1.4* betalogRR / sum(abs(betalogRR))
  }

  RR <- exp(V %*% betalogRR)
 
  Q <- (1-A)*Q0 + A*Q0*RR
  Y <- rpois(n, lambda = Q) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q, pY0 = Q0, pY1 =Q0*RR,  RR= RR )
  return(list(descr = "Data simulated from parametric linear model with known relative risk", beta_logRR = betalogRR, data = data, W = W, A = A, Y= Y, beta_A  = beta_A, beta_Y0 = beta_Y, link = "logistic"))
}
#' Simulate a dataset with known constant OR
#' @export
#' @param n Sample size
#' @param p Dimension of W
#' @param prop_active Proportion of active (nonzero coef) variables
sim.OR <- function(n=1500, p=2, prop_active = 1, formula_estimand =~1, formula_A = ~., formula_Y0W = ~., beta = NULL, beta_A = NULL, beta_Y = NULL) {
  W <- as.matrix(replicate(p, runif(n, -1,1)))
  colnames(W) <- paste0("W", 1:p)
  XA <- model.matrix(formula_A,as.data.frame(W))
  XY <- model.matrix(formula_Y0W,as.data.frame(W))
  V <- model.matrix(formula_estimand,as.data.frame(W))
  pA<- ncol(XA)
  pY <- ncol(XY)
  if(is.null(beta_A)) {
    activeA <- rbinom(pA, size = 1, prob =prop_active)
    beta_A <-  runif(pA, min=-1,max=1)
    beta_A <-   beta_A * activeA / sum(abs(beta_A*activeA))
    names(beta_A) <- colnames(XA)
  }
  if(is.null(beta_Y)) {
    activeY <- rbinom(pY, size = 1, prob =prop_active)
    beta_Y <- runif(pY, min=-1,max=1)
    beta_Y <- 0.7*beta_Y*activeY / sum(abs(beta_Y*activeY))
    names(beta_Y) <- colnames(XY)
  }
  g1 <- plogis( XA %*% beta_A)
  A <- rbinom(n, size = 1, prob = g1)
   
  
   
  if(!is.null(beta)) {
    betalogOR <- beta
  } else {
    if(ncol(V) > 1){
    betalogOR <- runif(ncol(V), min=-1,max=1)
    betalogOR <-   1.5*betalogOR / sum(abs(betalogOR))
    betalogOR <- c( betalogOR)
    } else{
      betalogOR <- 1
    }
  }
  
  logOR <- V%*%betalogOR 
  Q <- plogis(  A*logOR  +  XY %*% beta_Y)
  Q1 <- plogis(   1*logOR  +  XY %*% beta_Y)
  Q0 <- plogis(  XY %*% beta_Y)
  Y <- rbinom(n, size = 1, prob = Q) 
  
  data <- data.frame(W, A=A, Y=Y, pA1 = g1, pY = Q , pY1 = Q1, pY0 = Q0, OR = exp(logOR))
  return(list(descr = "Data simulated from parametric linear model with known odds ratio", beta_logOR = betalogOR, data = data, W = W, A = A, Y= Y, beta_A  = beta_A, beta_Y0 = beta_Y, link = "logistic"))
  
}


 