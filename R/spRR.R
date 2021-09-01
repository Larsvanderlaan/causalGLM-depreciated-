


#' Semiparametric targeted conditional relative risk/treatment effect using machine-learning
#' RR(W) := E[Y|A=1/W] / E[Y|A=0,W]
#' @param formula_logRR R-formula object specifying model for log relative risk
#' @param family_RR A R-family object specifying the link function for the log relative risk (gaussian/identity implies formula_logRR is directly modelling the log RR)
#' @param W A matrix of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y A nonnegative outcome variable. Can be binary, a count, or a continuous nonnegative variable.
#' @param sl3_Lrnr_A An optional sl3-Learner object to estimate P(A=1|W)
#' @param sl3_Lrnr_Y An optional sl3-Learner object to estimate nuisance conditional means E[Y|A=0,W] and E[Y|A=1,W]
#' @param weights A vector of optional weights.
#' @param smoothness_order Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param num_knots Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param max_degree Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#' @param fit_control Specification for default HAL learner (used if sl3 Learners not given). See spOR for use.
#'
spRR <- function(formula_logRR =  ~1, W, A, Y, family_RR = gaussian(),  sl3_Lrnr_A = NULL, sl3_Lrnr_Y = NULL, weights = NULL,  smoothness_order = 1, max_degree = 2, num_knots = c(15,5), fit_control = list()){
  fit_separate <- !is.null(sl3_Lrnr_Y) || family_RR$family != "gaussian" || family_RR$link != "identity"
  default_learner <- Lrnr_hal9001$new(smoothness_orders = smoothness_order, num_knots = num_knots, max_degree = max_degree, fit_control = fit_control )

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
  V <- model.matrix(formula_logRR , data = dat)

  # Estimate g
  data_A <- data.frame(W, A = A, weights = weights)
  task_A <- sl3_Task$new(data_A, covariates = colnames(W), outcome = "A", weights= "weights")
  sl3_Lrnr_A <- sl3_Lrnr_A$train(task_A)
  g1 <- sl3_Lrnr_A$predict(task_A)
  g1 <- bound(g1, 0.005)
  g0 <- 1- g1


  # Estimate part lin Q



  data_Y <- data.frame(W, A = A, Y=Y, weights = weights)
  task_Y <- sl3_Task$new(data_Y, covariates = c(colnames(W), "A"), outcome = "Y" , weights= "weights")
  data_Y1 <- data.frame(W, A = 1, Y=Y, weights = weights)
  task_Y1 <- sl3_Task$new(data_Y1, covariates = c(colnames(W), "A"), outcome = "Y", weights= "weights")
  data_Y0 <- data.frame(W, A = 0, Y=Y, weights = weights)
  task_Y0 <- sl3_Task$new(data_Y0, covariates = c(colnames(W), "A"), outcome = "Y", weights= "weights")

  sl3_Lrnr_Y <- sl3_Lrnr_Y$train(task_Y)
  Q <-  pmax(sl3_Lrnr_Y$predict(task_Y),5*1e-3)
  Q1 <-  pmax(sl3_Lrnr_Y$predict(task_Y1),5*1e-3)
  Q0 <-  pmax(sl3_Lrnr_Y$predict(task_Y0),5*1e-3)



  if(family_RR$family == "gaussian" & family_RR$link == "identity") {
    beta <- suppressWarnings(coef(glm.fit(V, Q1/Q0, family = poisson(), intercept = F)))
  } else {
    beta <- coef(glm.fit(V, log(Q1/Q0), family = family_RR, intercept = F))
  }
  link <- V %*% beta
  logRR <- as.vector(family_RR$linkinv(link))
  RR <- as.vector(exp(logRR))
  Q0 <- as.vector(as.vector(Q0))
  Q1 <- as.vector(Q0 * RR)
  Q <- as.vector(ifelse(A==1, Q1, Q0))
  print(range(Q1))
  print(range(Q0))
  print(range(RR))






  for(i in 1:100) {
    gradM <- family_RR$mu.eta(V%*%beta)*V
    mstar <- RR + (1-A)*1
    num <- gradM * ( RR * g1)
    denom <- RR * g1 + g0
    hstar <- - num/denom
    H <- (A*gradM  + hstar)

    EIF <- weights *   as.matrix(H * (Y-Q))

    linpred <- family_RR$linkfun(log(Q1/Q0))
    risk_function <- function(beta) {
      logQeps <- A*family_RR$linkinv(linpred + V%*%beta ) + log(Q0)+ hstar%*%beta
      loss <- exp(logQeps) - Y * logQeps
      loss <- weights*loss
      mean(loss)
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
    Q0 <- exp(log(Q0) + hstar %*% eps)
    RR <- exp(family_RR$linkinv(linpred +  V %*% eps))
    Q1 <- RR*Q0
    Q <- ifelse(A==1, Q1, Q0)
    beta <- coef(glm.fit(V, logRR, family = family_RR, intercept = F))

  }

  est <- beta
  se <- sqrt(diag(var(EIF)))
  Zvalue <- abs(sqrt(n) * est/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)

  ci <- cbind(est - 1.96*se/sqrt(n),est +1.96*se/sqrt(n) )
  out <- cbind(est, se, se/sqrt(n), Zvalue, ci, pvalue)
  colnames(out) <- c("coef", "se", "se/sqrt(n)", "Z-score", "CI_left", "CI_right", "p-value")
  out
}
