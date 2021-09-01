


#' Semiparametric targeted conditional average treatment effect estimation.
#' CATE(W) := E[Y|A=1,W] - E[Y|A=0,W]
#' @param formula_CATE R-formula object specifying model for CATE
#' @param family_CATE A R-family object specifying the link function for the CATE
#' @param W A matrix of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y An outcome variable (continuous or binary)
#' @param sl3_Lrnr_A An optional sl3-Learner object to estimate P(A=1|W)
#' @param sl3_Lrnr_Y An optional sl3-Learner object to estimate nuisance conditional means E[Y|A=0,W] and E[Y|A=1,W]
#' @param weights A vector of optional weights.
#' @param smoothness_order_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param num_knots_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param max_degree_Y0W Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param fit_control Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#'
#' @export
spCATE <- function(formula_CATE =  ~1, W, A, Y, family_CATE = gaussian(),  sl3_Lrnr_A = NULL, sl3_Lrnr_Y = NULL, weights = NULL,  smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(15,5), fit_control = list()){
  fit_separate <- !is.null(sl3_Lrnr_Y) || family_CATE$family != "gaussian" || family_CATE$link != "identity"
  default_learner <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, num_knots = num_knots_Y0W, max_degree = max_degree_Y0W, fit_control = fit_control )

  W <- as.matrix(W)
  A <- as.vector(A)
  Y <- as.vector(Y)
  n <- nrow(W)

  if(is.null(weights)) {
    weights <- rep(1,n)
  }
  fit_control$weights <- weights
  if(is.null(sl3_Lrnr_Y)) {
    sl3_Lrnr_Y <- default_learner
  }
  if(is.null(sl3_Lrnr_A)) {
    sl3_Lrnr_A <- default_learner
  }
   





  dat <-  as.data.frame(W)

  V <- model.matrix(formula_CATE , data = dat)

  # Estimate g
  data_A <- data.frame(W, A = A, weights = weights)
  task_A <- sl3_Task$new(data_A, covariates = colnames(W), outcome = "A", weights= "weights")
  sl3_Lrnr_A <- sl3_Lrnr_A$train(task_A)
  g1 <- sl3_Lrnr_A$predict(task_A)
  g0 <- 1- g1



  fit_separate <- !is.null(sl3_Lrnr_Y)

  # Estimate part lin Q
  if(!fit_separate){
    fit_Y <- fit_hal(X = as.matrix(W), X_unpenalized = as.matrix(A*V), Y = as.vector(Y), family = family_CATE, fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W)
    Q <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (A*V))
    Q0 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (0*V))
    Q1 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (1*V))
  } else {

    data_Y <- data.frame(W, A = A, Y=Y, weights = weights)
    task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(W), "A"), outcome = "Y" , weights= "weights")
    data_Y1 <- data.frame(W, A = 1, Y=Y, weights = weights)
    task_Y1 <- sl3_Task$new(data_Y1, covariates = c(colnames(W), "A"), outcome = "Y", weights= "weights")
    data_Y0 <- data.frame(W, A = 0, Y=Y, weights = weights)
    task_Y0 <- sl3_Task$new(data_Y0, covariates = c(colnames(W), "A"), outcome = "Y", weights= "weights")

    sl3_Lrnr_Y <- sl3_Lrnr_Y$train(task_Y)
    Q <-  sl3_Lrnr_Y$predict(task_Y)
    Q1 <-  sl3_Lrnr_Y$predict(task_Y1)
    Q0 <-  sl3_Lrnr_Y$predict(task_Y0)

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
    sigma2 <- Q*(1-Q)
    sigma21 <- Q1*(1-Q1)
    sigma20 <- Q0*(1-Q0)
  } else {
    X <- cbind(W,A)
    X0 <- cbind(W,rep(0,n))
    X1 <- cbind(W,rep(1,n))
    fit_Y <- fit_hal(X = X, , Y = (Y - Q)^2, family = "poisson", fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W)
    sigma2 <- predict(fit_Y, new_data =X)
    sigma20 <- predict(fit_Y, new_data = X0)
    sigma21 <- predict(fit_Y, new_data = X1)

    # data_sigma <- data.frame(W, A = A, Y=(Y-Q)^2, weights = weights)
    # task_sigma <- sl3_Task$new(data_sigma, covariates = c(colnames(W), "A"), outcome = "Y", weights = "weights")
    # data_sigma1 <- data.frame(W, A = 1, Y=(Y-Q)^2, weights = weights)
    # task_sigma1 <- sl3_Task$new(data_sigma1, covariates = c(colnames(W), "A"), outcome = "Y", weights = "weights")
    # data_sigma0 <- data.frame(W, A = 0, Y=(Y-Q)^2, weights = weights)
    # task_sigma0 <- sl3_Task$new(data_sigma0, covariates = c(colnames(W), "A"), outcome = "Y", weights = "weights")
    #
    # sl3_Lrnr_sigma <- sl3_Lrnr_sigma$train(task_sigma)
    # sigma2 <- sl3_Lrnr_sigma$predict(task_sigma)
    # sigma21 <- sl3_Lrnr_sigma$predict(task_sigma1)
    # sigma20 <- sl3_Lrnr_sigma$predict(task_sigma0)
    #sigma2 <- EY2 - Q^2
    #sigma20 <- EY2_0 - Q0^2
    #sigma21 <- EY2_1 - Q1^2
  }



  for(i in 1:100) {
    gradM <- family_CATE$mu.eta(V%*%beta)*V

    num <- gradM * ( g1/sigma21)
    denom <- (g0/ sigma20 + g1/sigma21)
    hstar <- - num/denom
    H <- (A*gradM  + hstar) /sigma2
    EIF <- weights * as.matrix(H * (Y-Q))

    linpred <- family_CATE$linkfun(Q1-Q0)

    risk_function <- function(beta) {
      loss <- weights*(Y - family_CATE$linkinv(A*linpred +    A*V %*% beta) - Q0 - hstar %*% beta)^2 / sigma2
      mean(loss)/2
    }
    suppressWarnings(hessian <-  optim(rep(0, ncol(V)),   fn = risk_function, hessian = T)$hessian)
    scale <- hessian

    #print(as.data.frame(hessian))

    #scale <- as.matrix(apply(gradM, 2, function(v) {colMeans_safe(weights*(A*gradM  + hstar) *  A*gradM * v /sigma2  )}) )
    #print(as.data.frame(scale))
    #stop("d")
    scaleinv <- solve(scale)
    EIF <-  EIF %*%   scaleinv


    scores <- colMeans(EIF)
    direction_beta <- scores/sqrt(mean(scores^2))
    print(max(abs(scores)))
    if(max(abs(scores)) <= 1/n) {
      break
    }
    linpred <- family_CATE$linkfun(Q1-Q0)
    risk_function <- function(eps) {

      loss <- weights*(Y - family_CATE$linkinv(A*linpred +  eps * A*V %*%direction_beta) - Q0 - eps*hstar %*% direction_beta)^2 / sigma2
      mean(loss)
    }

    optim_fit <- optim(
      par = list(epsilon = 0.01), fn = risk_function,
      lower = 0, upper = 0.01,
      method = "Brent"
    )
    eps <-  direction_beta * optim_fit$par
    Q0 <- as.vector(Q0 + hstar %*% eps)
    CATE <- family_CATE$linkinv(linpred +  V %*% eps)
    beta <- coef(glm.fit(V, CATE, family = family_CATE, intercept = F))
    link <- as.vector(V %*% beta)
    CATE <- family_CATE$linkinv(link)
    Q <- as.vector(A*CATE + Q0)
    Q1 <- as.vector(CATE + Q0)
  }

  est <- beta
  se <- sqrt(diag(var(EIF)))
  Zvalue <- abs(sqrt(n) * est/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)

  ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
  out <- cbind(est, se/sqrt(n), se,    ci,  Zvalue,pvalue)
  colnames(out) <- c("coef", "se/sqrt(n)", "se", "CI_left", "CI_right",  "Z-score", "p-value")
  output <- (list(coefs = out, var_mat = var(EIF)))
  class(output) <- c("spCATE", "causalGLM")
  output
}
