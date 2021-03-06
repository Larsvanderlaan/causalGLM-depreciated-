


#' Semiparametric targeted conditional relative risk/treatment effect using machine-learning
#' RR(W) := E[Y|A=1/W] / E[Y|A=0,W]
#' @param formula_logRR R-formula object specifying model for log relative risk
#' @param family_RR A R-family object specifying the link function for the log relative risk (gaussian/identity implies formula_logRR is directly modelling the log RR)
#' @param W A matrix of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y A nonnegative outcome variable. Can be binary, a count, or a continuous nonnegative variable.
#' @param pool_A_when_training Whether to estimate E[Y|A=1,W] and E[Y|A=0,W] separately or pool the regression.
#' If TRUE, the design matrix passed to the regression algorithm/learner is cbind(W,A*V) where V is the design matrix specified by the argument \code{formula_logRR}.
#' 
#' @param sl3_Learner_A An optional sl3-Learner object to estimate P(A=1|W)
#' @param sl3_Learner_Y An optional sl3-Learner object to estimate nuisance conditional means E[Y|A=0,W] and E[Y|A=1,W] (either pooled or separately depending on \code{pool_A_when_training})
#' If pool_A_when_training = T, it is recommended to specify the sl3 internal family object of the learner to equal `poisson()`.`
#' @param weights A vector of optional weights.
#' @param smoothness_order_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param num_knots_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param max_degree_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param fit_control Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @export
spRR <- function(formula_logRR =  ~1, W, A, Y, family_RR = gaussian(), pool_A_when_training = TRUE, full_fit_as_offset = TRUE, sl3_Learner_A = NULL, sl3_Learner_Y = NULL, weights = NULL,  smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(15,5), fit_control = list(), return_competitor = F){
  fit_separate <- !is.null(sl3_Learner_Y) || family_RR$family != "gaussian" || family_RR$link != "identity"
  default_learner <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, num_knots = num_knots_Y0W, max_degree = max_degree_Y0W, fit_control = fit_control )
  
  W <- as.matrix(W)
  A <- as.vector(A)
  Y <- as.vector(Y)
  n <- nrow(W)
  
  if(is.null(weights)) {
    weights <- rep(1,n)
  }
  fit_control$weights <- weights
   
  if(is.null(sl3_Learner_A)) {
    sl3_Learner_A <- default_learner
  }
  
  
  
  
  
  
  dat <-  as.data.frame(W)
  V <- model.matrix(formula_logRR , data = dat)
  
  # Estimate g
  data_A <- data.frame(W, A = A, weights = weights)
  task_A <- sl3_Task$new(data_A, covariates = colnames(W), outcome = "A", weights= "weights")
  sl3_Learner_A <- sl3_Learner_A$train(task_A)
  g1 <- sl3_Learner_A$predict(task_A)
  g1 <- as.vector(bound(g1, 0.005))
  g0 <- 1- g1
  
  
  # Estimate part lin Q
  binary <- all(Y %in% c(0,1))
  
  if(binary){
    outcome_type = "binomial"
  } else {
    outcome_type = "continuous"
  }
  if(is.null(sl3_Learner_Y)) {
     
    fit_control$weights <- weights
    fit_Y <- fit_hal(X = as.matrix(W), X_unpenalized = as.matrix(A*V), Y = as.vector(Y), family = "poisson", return_x_basis = T, fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W)
    Q <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (A*V))
    Q0 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (0*V))
    Q1 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (1*V))
  } else {
    if(pool_A_when_training) {
      Vtmp <- V
      colnames(Vtmp) <- paste0("V", 1:ncol(V))
      X <- cbind(W, A*Vtmp)
      X1 <- cbind(W, Vtmp)
      X0 <- cbind(W, 0*Vtmp)
      data_Y <- data.table(X, A = A, Y=Y, weights = weights)
      
     
      task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(X)), outcome = "Y" , weights= "weights", outcome_type = outcome_type)
      data_Y1 <- data.table(X1, A = 1, Y=Y, weights = weights)
      task_Y1 <- sl3_Task$new(data_Y1, covariates = c(colnames(X)), outcome = "Y", weights= "weights")
      data_Y0 <- data.table(X0, A = 0, Y=Y, weights = weights)
      task_Y0 <- sl3_Task$new(data_Y0, covariates = c(colnames(X)), outcome = "Y", weights= "weights")
        
      sl3_Learner_Y <- sl3_Learner_Y$train(task_Y)
      Q <-  pmax(sl3_Learner_Y$predict(task_Y),0.01)
      Q1 <-  pmax(sl3_Learner_Y$predict(task_Y1),0.01)
      Q0 <-  pmax(sl3_Learner_Y$predict(task_Y0),0.01)
      
    } else {
      data_Y <- data.table(W, Y=Y, weights = weights)
      task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(W)), outcome = "Y" , weights= "weights", outcome_type = outcome_type)
      lrnr_Y0 <- sl3_Learner_Y$train(task_Y)$train(task_Y[A==0])
      lrnr_Y1 <- sl3_Learner_Y$train(task_Y)$train(task_Y[A==1])
      Q1 <-  pmax(lrnr_Y1$predict(task_Y),0.01)
      Q0 <-  pmax(lrnr_Y0$predict(task_Y),0.01)
      Q <- ifelse(A==1, Q1, Q0)
    }
  }
  beta <- suppressWarnings(coef(glm.fit(V, log(Q1/Q0), family = gaussian(), intercept = F)))
  
   
  if(full_fit_as_offset ) {
    beta <- coef(glm.fit(A*V, Y,  offset = log(Q), family = poisson(), intercept = F))
    Q <- exp(log(Q) + A*V %*%beta)
    Q1 <- exp(log(Q1) + V %*%beta)
    
  }
  beta <- suppressWarnings(coef(glm.fit(V, log(Q1/Q0), family = gaussian(), intercept = F)))
  
  Q0 <- as.vector(pmax(Q0,0.01))
  Q1 <- as.vector(pmax(Q1,0.01))
  Q <- as.vector(ifelse(A==1, Q1, Q0))
  
   
  if(family_RR$family == "gaussian" & family_RR$link == "identity") {
    beta <- suppressWarnings(coef(glm.fit(V, log(Q1/Q0), family = gaussian(), intercept = F)))
  } else {
    beta <- coef(glm.fit(V, log(Q1/Q0), family = family_RR, intercept = F))
  }
  link <- V %*% beta
  logRR <- as.vector(family_RR$linkinv(link))
  RR <- as.vector(exp(logRR))
  Q0 <- pmax(as.vector(as.vector(Q0)),0.001)
  Q1 <- pmax(as.vector(Q0 * RR),0.001)
  Q <- as.vector(ifelse(A==1, Q1, Q0))
  
  
  
  ###  Est EQN
  
  risk_function <- function(beta) {
    gradM <- family_RR$mu.eta(V%*%beta)*V
    RR <- as.vector(exp(V %*% beta))
    mstar <- RR + (1-A)*1
    num <- gradM * ( RR * g1)
    denom <- RR * g1 + g0
    hstar <- - num/denom
    H <- (A*gradM  + hstar)
    Q <- pmax(A*RR*Q0 + (1-A)*Q0, 0.001)
    EIF <- weights *   as.matrix(H * (Y-Q))
    
    sds <- apply(EIF,2,sd)
    sds <- 1/sds
    sds <- sds/sum(sds)
    sds <- 1
    (sum(sds*(colMeans(EIF)^2)))
  }
  (one_step <-  optim(rep(0, ncol(V)),   fn = risk_function,  method = "BFGS" ))
   
  
  one_step <- one_step$par
 
  
 
   EIF_init <- NULL
  for(i in 1:200) {
    gradM <- family_RR$mu.eta(V%*%beta)*V
    
    mstar <- RR + (1-A)*1
    num <- gradM * ( RR * g1)
    denom <- RR * g1 + g0
    hstar <- - num/denom
    H <- (A*gradM  + hstar)
    
    EIF <- weights *   as.matrix(H * (Y-Q))
    
    
    scale <- apply(V,2, function(v) {
      colMeans_safe(weights*V*v*g1*g0*RR/(g1*RR + g0)^2 *(Y-Q) + H*(A*v*Q))
    })
    
    
    linpred <- family_RR$linkfun(log(Q1/Q0))
    risk_function <- function(beta) {
      logQeps <- A*family_RR$linkinv(linpred + V%*%beta ) + log(Q0)+ hstar%*%beta
      loss <- exp(logQeps) - Y * logQeps
      loss <- weights*loss
      mean(loss)
    }
    
    if(family_RR$family == "gaussian" && family_RR$link == "identity") {
      scale <- apply(V,2, function(v) {
        colMeans_safe(weights*V*v*g1*g0*RR/(g1*RR + g0)^2 *(Y-Q) + H*(A*v*Q))
      })
       
    } else {
      suppressWarnings(scale <-  optim(rep(0, ncol(V)),   fn = risk_function, hessian = T)$hessian)
    }
    #suppressWarnings(scale <-  optim(rep(0, ncol(V)),   fn = risk_function, hessian = T)$hessian)
    
     
    
    
    #print(as.data.frame(hessian))
    
    #scale <- as.matrix(apply(gradM, 2, function(v) {colMeans_safe(weights*(A*gradM  + hstar) *  A*gradM * v /sigma2  )}) )
    #print(as.data.frame(scale))
    #stop("d")
    scaleinv <- solve(scale)
    
    EIF <-  EIF %*%   scaleinv
    if(is.null(EIF_init)){
      EIF_init <- EIF
    }
    
    scores <- colMeans(EIF)
    print(scores)
    direction_beta <- scores/sqrt(mean(scores^2))
    #print("scores")
    #print(scores)
     
    EIF_weights <- 1/(diag(var(EIF_init)))
    EIF_weights <- pmin(EIF_weights, (sqrt(n)/log(n))^2)
    direction_beta <- scores * sqrt(EIF_weights) 
    direction_beta <- direction_beta/sqrt(mean(direction_beta^2)) 
    scores <- sqrt(sum(scores^2*EIF_weights))
    
    
    if(abs(scores) <= 0.5/sqrt(n)/log(n)){
     
      break
    }
    
    linpred <- family_RR$linkfun(log(Q1/Q0))
    risk_function <- function(eps) {
      logQeps <- A*family_RR$linkinv(linpred + eps * V%*%direction_beta ) + log(Q0)+ eps * hstar%*%direction_beta
      loss <- exp(logQeps) - Y * logQeps
      loss <- weights*loss
      mean(loss)
    }
    
    
    optim_fit <- optim(
      par = list(epsilon = 0.01), fn = risk_function,
      lower = 0, upper = 0.01,
      method = "Brent"
    )
    eps <-  direction_beta * optim_fit$par
    Q0 <-  as.vector(pmax(exp(log(Q0) + hstar %*% eps), 0.005))
    RR <- as.vector(exp(family_RR$linkinv(linpred +  V %*% eps)))
    logRR <- log(RR)
    Q1 <-  RR*Q0 
    Q <- ifelse(A==1, Q1, Q0)
     beta <- coef(glm.fit(V, logRR, family = family_RR, intercept = F))
    # logRR <- as.vector(V%*% beta)
    # RR <- exp(logRR)
    # Q1 <- pmax(RR*Q0, 0.001)
    # Q <- ifelse(A==1, Q1, Q0)
  }
   
  beta <- coef(glm.fit(V, logRR, family = family_RR, intercept = F))
  
  est <- beta
  var_mat <- var(EIF_init)
  
  RR_linkinv <- function(x) {exp(family_RR$linkinv(x))}
  
  se <- sqrt(diag(var_mat ))
  Zvalue <- abs(sqrt(n) * est/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)
  
  ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
  out <- cbind(est, se/sqrt(n), se,    ci,  Zvalue,pvalue)
  colnames(out) <- c("coefs", "se/sqrt(n)", "se", "CI_left", "CI_right",  "Z-score", "p-value")
  
  output <- list(coefs = out, var_mat = var(EIF), n=n,  formula = formula_logRR, linkinv = RR_linkinv, link_type = "Maps linear predictor to RR")
  class(output) <- c("spRR", "causalGLM")
  output
  
  if(return_competitor){
  est <- one_step
  se <- sqrt(diag(var(EIF_init)))
  Zvalue <- abs(sqrt(n) * est/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)
  ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
  out <- cbind(est, se/sqrt(n), se,    ci,  Zvalue,pvalue)
  colnames(out) <- c("coefs", "se/sqrt(n)", "se", "CI_left", "CI_right",  "Z-score", "p-value")
  output$coefs1 <- out
  }
  
  
  output
}
