


#' Semiparametric targeted conditional average treatment effect estimation.
#' CATE(W) := E[Y|A=1,W] - E[Y|A=0,W]
#' @param formula_CATE R-formula object specifying model for CATE
#' @param family_CATE A R-family object specifying the link function for the CATE
#' @param W A matrix of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y An outcome variable (continuous or binary)
#' @param sl3_Learner_A An optional sl3-Learner object to estimate P(A=1|W)
#' @param sl3_Learner_Y An optional sl3-Learner object to estimate nuisance conditional means E[Y|A=0,W] and E[Y|A=1,W]
#' @param weights A vector of optional weights.
#' @param smoothness_order_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param num_knots_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param max_degree_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param fit_control Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#'
#' @export
spCATE <- function(formula_CATE =  ~1, W, A, Y, family_CATE = gaussian(), pool_A_when_training = T, constant_variance = FALSE,  return_ATE = TRUE, sl3_Learner_A = NULL, sl3_Learner_Y = NULL, weights = NULL,  smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(15,5), max_degree_sigma =1, num_knots_sigma = ifelse(ncol(W)>=25, 1, 10), fit_control = list(), return_competitor = FALSE){
   
  fit_separate <- !pool_A_when_training  
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

  V <- model.matrix(formula_CATE , data = dat)

  # Estimate g
  data_A <- data.frame(W, A = A, weights = weights)
  task_A <- sl3_Task$new(data_A, covariates = colnames(W), outcome = "A", weights= "weights")
  sl3_Learner_A <- sl3_Learner_A$train(task_A)
  g1 <- sl3_Learner_A$predict(task_A)
  g0 <- 1- g1


 
  # Estimate part lin Q
  if(is.null(sl3_Learner_Y) & !(family_CATE$family == "gaussian" && family_CATE$link == "identity")) {
    sl3_Learner_Y <- default_learner
  }
  if(is.null(sl3_Learner_Y)){
     
    
     
    fit_Y <- fit_hal(X = as.matrix(W), X_unpenalized = as.matrix(A*V), Y = as.vector(Y), family = "gaussian", fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W)
    Q <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (A*V))
    Q0 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (0*V))
    Q1 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (1*V))
  } else {
    X <- W
    X1 <- W
    X0 <- W
    data_Y <- data.frame(X, A = A, Y=Y, weights = weights)
    task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(X), "A"), outcome = "Y" , weights= "weights")
    task_Y0 <- sl3_Task$new(data_Y, covariates =  colnames(X) , outcome = "Y" , weights= "weights") 
   
    if(fit_separate && family_CATE$family == "gaussian" && family_CATE$link == "identity") {
     
      sl3_Learner_Y <- sl3_Learner_Y$train(task_Y0[A==0])
      Q0 <-  sl3_Learner_Y$predict(task_Y0)
      beta <- coef(glm.fit(A*V, Y, family = gaussian(), offset = Q0, intercept = F))
      Q1 <- Q0 + V %*% beta
      Q <- ifelse(A==1, Q1, Q0)
    } else if(fit_separate) {
      sl3_Learner_Y0 <- sl3_Learner_Y$train(task_Y0[A==0])
      sl3_Learner_Y1 <- sl3_Learner_Y$train(task_Y0[A==1] )
      Q1 <-  sl3_Learner_Y1$predict(task_Y0)
      Q0 <-  sl3_Learner_Y0$predict(task_Y0)
      Q <- ifelse(A==1, Q1,Q0)
       
    } else {
      
      Vtmp <- V
      colnames(Vtmp) <- paste0("V", 1:ncol(V))
      X <- cbind(W, A*Vtmp)
      X1 <- cbind(W, Vtmp)
      X0 <- cbind(W, 0*Vtmp)
      
      data_Y <- data.table(X,   Y=Y, weights = weights)
      task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(X) ), outcome = "Y" , weights= "weights")
      data_Y1 <- data.table(X1,    Y=Y, weights = weights)
      task_Y1 <- sl3_Task$new(data_Y1, covariates = c(colnames(X) ), outcome = "Y", weights= "weights")
      data_Y0 <- data.table(X0,  Y=Y, weights = weights)
      task_Y0 <- sl3_Task$new(data_Y0, covariates = c(colnames(X) ), outcome = "Y", weights= "weights")
      sl3_Learner_Y <- sl3_Learner_Y$train(task_Y )
      Q   <-  sl3_Learner_Y$predict(task_Y)
      Q1 <-  sl3_Learner_Y$predict(task_Y1)
      Q0 <-  sl3_Learner_Y$predict(task_Y0)
    }
     
  }

  
  beta <- coef(glm.fit(V, Q1-Q0, family = family_CATE, intercept = F))
  link <- V %*% beta
  CATE <- family_CATE$linkinv(link)
  Q0 <- as.vector(Q0)
  Q <- as.vector(A*CATE + Q0)
  Q1 <- as.vector(CATE + Q0)




  # Estimate var
   
 
  binary <- all(Y %in% c(0,1))
  if(binary) {
    Qtmp <- bound(Q,0.005)
    Qtmp1 <- bound(Q1,0.005)
    Qtmp0 <- bound(Q0,0.005)
    sigma2 <- Qtmp*(1-Qtmp)
    sigma21 <- Qtmp1*(1-Qtmp1)
    sigma20 <- Qtmp0*(1-Qtmp0)
  } else {
    if(constant_variance) {
      sigma2 <- mean((Y - Q)^2)
      sigma2 <- rep(sigma2, n)
      sigma21 <- mean((Y - Q1)^2)
      sigma21 <- sigma2
      sigma20 <- mean((Y - Q0)^2)
      sigma20 <- sigma2
      
    } else {
      X <- cbind(W,A, A*W)
      X0 <- cbind(W,rep(0,n), 0*W)
      X1 <- cbind(W,rep(1,n), W)
      fit_Y <- fit_hal(X = X, , Y = (Y - Q)^2, family = "poisson", fit_control = fit_control, smoothness_orders = 1, max_degree = max_degree_sigma, num_knots = num_knots_sigma)
      sigma2 <- predict(fit_Y, new_data =X)
      sigma20 <- predict(fit_Y, new_data = X0)
      sigma21 <- predict(fit_Y, new_data = X1)
       
    }
    
  }

   
  #one_step <- mean((A-g1)*(Y-Q0))/ mean((A-g1)*A)
  
  
  
  
  risk_function <- function(beta) {
    gradM <- family_CATE$mu.eta(V%*%beta)*V
    num <- gradM * ( g1/sigma21)
    denom <- (g0/ sigma20 + g1/sigma21)
    hstar <- - num/denom
    H <- (A*gradM  + hstar) /sigma2
    EIF <- weights * as.matrix(H * (Y -  A*as.vector(V%*%beta) - Q0))
    sds <- apply(EIF,2,sd)
    sds <- 1/sds
    sds <- 1
    
    
    (sum(sds*(colMeans(EIF)^2)))
  }
   (one_step <-  optim(rep(0, ncol(V)),   fn = risk_function, method = "BFGS"))
 
  one_step <- one_step$par
   
  
  
  
  
  
  hessian <- NULL
  EIF_init <- NULL
  for(i in 1:10) {
    gradM <- family_CATE$mu.eta(V%*%beta)*V

    num <- gradM * ( g1/sigma21)
    denom <- (g0/ sigma20 + g1/sigma21)
    hstar <- - num/denom
     
    H <- (A*gradM  + hstar) /sigma2
    EIF <- weights * as.matrix(H * (Y-Q))
     
    linpred <- family_CATE$linkfun(Q1-Q0)
    if(T || is.null(hessian)){
    risk_function <- function(beta) {
      loss <- weights*(Y - family_CATE$linkinv(A*linpred +    A*V %*% beta) - Q0 - hstar %*% beta)^2 / sigma2
      mean(loss)/2
    }
    (hessian <-  optim(rep(0, ncol(V)),   fn = risk_function, hessian = T , method = "BFGS")$hessian)
    }
    scale <- hessian
    
    #print(as.data.frame(hessian))

    #scale <- as.matrix(apply(gradM, 2, function(v) {colMeans_safe(weights*(A*gradM  + hstar) *  A*gradM * v /sigma2  )}) )
    #print(as.data.frame(scale))
    #stop("d")
    scaleinv <- solve(scale)
    EIF <-  EIF %*%   scaleinv
    if(is.null(EIF_init)) {
      EIF_init <- EIF
    }

    scores <- colMeans(EIF)
    
    direction_beta <- scores/sqrt(mean(scores^2))
     #print(scores)
    if(max(abs(scores)) <= 1e-10) {
      break
    }
    linpred <- family_CATE$linkfun(Q1-Q0)
    risk_function <- function(eps) {
      #loss <- weights*(Y - family_CATE$linkinv(A*linpred + A*V %*%eps) - Q0 - hstar %*% eps)^2 / sigma2
      loss <- weights*(Y - family_CATE$linkinv(A*linpred +  eps * A*V %*%direction_beta) - Q0 - eps*hstar %*% direction_beta)^2 / sigma2
      mean(loss)
    }

    optim_fit <- optim(
      par = list(epsilon = 0)
        , fn = risk_function ,# method = "BFGS"
     lower = 0, upper = 1,
      method = "Brent"
    )
    eps <- direction_beta * optim_fit$par #optim_fit$par #direction_beta * optim_fit$par
    #print("eps")
    #print(eps)
    Q0 <- as.vector(Q0 + hstar %*% eps)
    if(binary) {
      Q0 <- bound(Q0,0)
    }
    CATE <- family_CATE$linkinv(linpred +  V %*% eps)
    beta <- coef(glm.fit(V, CATE, family = family_CATE, intercept = F))
    link <- as.vector(V %*% beta)
    CATE <- family_CATE$linkinv(link)
    Q <- as.vector(A*CATE + Q0)
    Q1 <- as.vector(CATE + Q0)
    
  }
  
   
  
  linkinv <- family_CATE$linkinv
  
  est <- beta
  se <- sqrt(diag(var(EIF_init)))
  se_keep <- se
  Zvalue <- abs(sqrt(n) * est/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)

  ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
  out <- cbind(est, se/sqrt(n), se,    ci,  Zvalue,pvalue)
  colnames(out) <- c("coefs", "se/sqrt(n)", "se", "CI_left", "CI_right",  "Z-score", "p-value")
  output <- list(coefs = out, ATE = NULL, var_mat = var(EIF_init), n=n, formula = formula_CATE, linkinv = linkinv, link_type = "Maps linear predictor to CATE", EIF = EIF)
  class(output) <- c("spCATE", "causalGLM")
  
  
  if(return_ATE & family_CATE$link == "identity") {
    
     
    
    EIF_ATE <- rowMeans(EIF_init%*%t(V))
    
    se <- sd(EIF_ATE)
    est <-  mean(V%*%beta)
 
    ci <- c(est - 1.96*se/sqrt(n), est + 1.96*se/sqrt(n))
    Zvalue <- abs(sqrt(n) * est/se)
    pvalue <- signif(2*(1-pnorm(Zvalue)),5)
    out <- matrix(c(est,se/sqrt(n) ,se, ci, Zvalue, pvalue), nrow=1)
    
    colnames(out) <- c("ATE",  "se/sqrt(n)","se", "CI_left", "CI_right", "Z-score", "p-value")
    output$ATE <- out
  }
  
  
  if(return_competitor) {
    est <- one_step
    se <- se_keep
    Zvalue <- abs(sqrt(n) * est/se)
    pvalue <- signif(2*(1-pnorm(Zvalue)),5)
    ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
    out <- cbind(est, se/sqrt(n), se,    ci,  Zvalue,pvalue)
    colnames(out) <- c("coefs", "se/sqrt(n)", "se", "CI_left", "CI_right",  "Z-score", "p-value")
    output$coefs1 <- out
   
  }
  output
   
}
