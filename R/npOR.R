

#' Nonparametric Targeted inference for the conditional odds ratio with outcome missingness.
#'
#' Nonparametric Targeted estimates and inference for the conditional odds ratio with a user-specified parametric working-model with informative outcome missingness due to `A` and `W`.
#' This version also allows for the outcome to be missing-at-random conditional on A,W.
#' Note this version is totally nonparametric and requires no assumptions on the functional form of the conditional odds ratio for correct inference.
#' The user-specified working model is used to find a best-approximation (relative to the working-model) of the true conditional odds ratio.
#' This function can be viewed as a nonparametric version of the partially-linear logistic regression model where one uses the working model `logit(P(Y=1|A,W)) = b* A f(W) +  logit(P(Y=1|A=0,W)` for a user-specified parametric function `f(W)` to approximate the true function.
#' Nonparametrically correct and efficient inference is given for this best approximation of the conditional odds ratio.
#'
#' @param working_formula An working-model R formula object describing the functional form of the approximation of the conditional log odds ratio as a fnction of `W`.
#' @param W A named matrix of baseline covariates
#' @param A A binary vector with values in (0,1) encoding the treatment assignment
#' @param Y A binary outcome variable with values in (0,1)
#' @param Delta A binary vector that takes the value 1 if `Y` is osberved/not-missing and 0 otherwise.
#' @param weights An optional vector of weights for each observation. Use with caution. This can lead to invalid inferences if included naively.
#' @param W_new An optional matrix of new values of the baseline covariates `W` at which to predict odds ratio.
#' @param glm_formula_Y0W (Not recommended). An optional R formula object describing the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`
#' @param sl3_learner_Y0W
#' @param glm_formula_OR (Not recommended). An optional R formula object describing the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`
#' @param sl3_learner_OR
#' @param glm_formula_A (Not recommended). An optional R formula object describing the functional form of P(A=1|W). If provided, \code{glm} is used for the fitting. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param  sl3_learner_A An optional \code{tlverse/sl3} learner object used to estimate P(A=1|W).
#'  If both \code{sl3_learner_A} and \code{glm_formula_A} are not provided, a default learner is used (Lrnr_hal9001).
#' @param glm_formula_Delta (Not recommended). An optional R formula object describing the functional form of P(Delta=1|A,W) to fit with glm. If provided, it is estimated using glm. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param sl3_learner_Delta An optional \code{tlverse/sl3} learner object used to estimate P(Delta=1|A,W).
#' @param sl3_learner_default A default sl3 Learner to be used if neither a glm formula or sl3 learner is provided for one of the nuisance functions.
#' By default, Lrnr_hal9001 is used.
#' @export
#'
npOR  <- function(working_formula = logOR~1, W, A, Y, Delta = NULL, weights = NULL, W_new = W,  glm_formula_Y0W = NULL, sl3_learner_Y0W = NULL,  glm_formula_OR = NULL, sl3_learner_OR = NULL, glm_formula_A = NULL, sl3_learner_A = NULL, glm_formula_Delta = NULL, sl3_learner_Delta = NULL, sl3_learner_default = Lrnr_hal9001_custom$new(max_degree =2, smoothness_orders = 1, num_knots = c(20,10)) ) {

  W <- as.matrix(W)
  if(is.null(Delta)) {
    Delta <- rep(1, nrow(W))
  }

  n <- nrow(W)
  formula <- as.formula(working_formula)
  V <- model.matrix(formula, data = as.data.frame(W))
  V_new <- model.matrix(formula, data = as.data.frame(W_new))
  if(is.null(weights)) {
    weights <- rep(1, nrow(W))
  }

  keep <- Delta==1

  ################################################################################################
  #### Learn P(Y=1|A=0,X)   #########################
  ################################################################################################
  if(is.null(sl3_learner_Y0W)) {
    sl3_learner_Y0W <- sl3_learner_default
  }
  if(is.null(glm_formula_Y0W)){
    # If no formula use HAL and respect model constraints.
    data <-  data.frame(W,Y)
    data$weights <- weights

    task_Y0W <- sl3_Task$new(data, covariates = colnames(W), outcome = "Y", weights = "weights" )
    sl3_learner_Y0W <- sl3_learner_Y0W$train(task_Y0W[keep & A==0])
    Q0 <- sl3_learner_Y0W$predict(task_Y0W)

  } else {
    # Use glm if formula supplied
    X <- model.matrix(glm_formula_Y0W, data = as.data.frame(W))

    fit_Y <- glm.fit(X,Y, weights = weights*keep*(A==0), family = binomial(), intercept = F)
    cfs <- coef(fit_Y)
    Q0 <- as.vector(plogis(X%*%cfs))
  }

  ################################################################################################
  #### Learn nonparametric odds ratio   #########################
  ################################################################################################
  if(is.null(sl3_learner_OR)) {
    sl3_learner_OR <- sl3_learner_default
  }
  if(is.null(glm_formula_OR)){
    # If no formula use HAL and respect model constraints.
    data <-  data.frame(W)
    covariates <- colnames(data)
    data$Y <- Y
    data$weights <- weights
    data$offset <- Q0

    task_Y1W <- sl3_Task$new(data, covariates = covariates, outcome = "Y" , weights = "weights", offset = "offset", outcome_type = "binomial")
    sl3_learner_OR <- sl3_learner_OR$train(task_Y1W[keep & A==1])
    Q1 <- sl3_learner_OR$predict(task_Y0W)

  } else {
    # Use glm if formula supplied
    X <- model.matrix(glm_formula_OR, data = as.data.frame(W))

    fit_Y <- glm.fit(X,Y, weights = weights*keep*(A==1), offset = qlogis(Q0), family = binomial(), intercept = F)
    cfs <- coef(fit_Y)
    Q1 <- as.vector(plogis(qlogis(Q0) + X%*%cfs))
  }

  Q <- ifelse(A==1, Q1, Q0)

  ################################################################################################
  #### Learn working-model odds ratio   #########################
  ################################################################################################
  fit_Y <- suppressWarnings(glm.fit(A*V,Q, offset = qlogis(Q0), weights = weights,family = binomial(), intercept = F))
  beta <- coef(fit_Y)
  logORbeta <- V%*%beta
  Q1beta <- as.vector(plogis(qlogis(Q0) + logORbeta ))
  ################################
  #### Learn P(A=1|X) ############
  ################################
  # Default sl3 learner
  # If no glm formula then use sl3 learner.
  if(is.null(sl3_learner_A)) {
    sl3_learner_A <- sl3_learner_default
  }
  if(is.null(glm_formula_A)) {

    data_A <- data.frame(W, A)
    covariates_A <- paste0("W", 1:ncol(W))
    colnames(data_A) <- c(covariates_A, "A")
    data_A$weights <- weights
    task_A <- sl3_Task$new(data_A, covariates = covariates_A, outcome = "A", outcome_type = "binomial", weights = "weights" )
    fit_A <- sl3_learner_A$train(task_A)
    g1 <- fit_A$predict(task_A)
  } else {
    # use glm is formula supplied
    W_g <- model.matrix(glm_formula_A, data = as.data.frame(W))
    fit_A <- glm.fit(W_g,A, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_A)
    g1 <- as.vector(as.vector(plogis(W_g%*%cfs)))
  }
  g0 <- 1-g1

  ################################
  #### Learn P(Delta=1|A,X) ############
  ################################
  if(is.null(sl3_learner_Delta)) {
    sl3_learner_Delta <- sl3_learner_default
  }
  if(is.null(glm_formula_Delta)) {

    data_Delta <- data.frame(W, A)
    covariates_Delta  <- c(paste0("W", 1:ncol(W)), "A")
    colnames(data_Delta ) <- covariates_Delta
    data_Delta$Delta <- Delta
    data_Delta$weights <- weights
    task_Delta <- sl3_Task$new(data_Delta, covariates = covariates_Delta, outcome = "Delta", outcome_type = "binomial", weights = "weights" )
    fit_Delta <- sl3_learner_Delta$train(task_Delta)
    G <- fit_Delta$predict(task_Delta)

    G <- pmax(G, 0.005)
  } else {
    # use glm is formula supplied
    W_np <- model.matrix(glm_formula_Delta, data = as.data.frame(cbind(W)))
    X <- cbind(W_np, A*W_np)
    X1 <- cbind(W_np, 1*W_np)
    X0 <- cbind(W_np, 0*W_np)
    fit_Delta<- glm.fit(X,Delta, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_Delta)
    G <- as.vector(plogis(X%*%cfs))
    G <- pmax(G, 0.005)
  }

  ################################
  ##### Targeting Step ###########
  ################################



  converged_flag <- FALSE
  for(i in 1:50) {
    sigma02 <- as.vector(Q0*(1-Q0))
    sigma1beta2 <- as.vector(Q1beta*(1-Q1beta))
    omega <- as.vector((g0*sigma02 + g1*sigma1beta2)/(g0*sigma02))

    h_star <- V* as.vector(-(g1*sigma1beta2) / (g1*sigma1beta2 + (1-g1)*sigma02))
    H_star <- (A*V + h_star)
    H_star1 <- (V + h_star)
    H_star0 <-  h_star
    offset <- qlogis(Q)
    scale <- apply(V,2, function(v){colMeans_safe(weights*(g1 * sigma1beta2 * v*V))})
    scale_inv <- solve(scale)



    EIFY <- weights*(Delta/G)*omega* (H_star %*% scale_inv) * (Y - Q)
    var_scaled <- as.matrix(var(EIFY))

    score <- sum(abs(colMeans_safe(EIFY) ))


    if(abs(score) <= min(0.5,mean(sqrt(diag(var_scaled))))/sqrt(n)/log(n)){
      converged_flag <- TRUE
      print("converged")
      break
    }
    offset <- qlogis(Q)
    eps <- coef(glm(Y~X-1, family = binomial(), weights = weights*(Delta/G)*omega, offset = offset, data = list(Y = Y, X = H_star)))
    Q <- as.vector(plogis(offset +  H_star %*% eps))
    Q0 <- as.vector(plogis(qlogis(Q0) +  H_star0 %*% eps ))
    Q1 <- as.vector(plogis(qlogis(Q1) +  H_star1 %*% eps))

    ################################################################################################
    #### Learn working-model odds ratio   #########################
    ################################################################################################
    fit_Y <- suppressWarnings(glm.fit(A*V,Q, offset = qlogis(Q0), weights = weights,family = binomial(), intercept = F))
    beta <- coef(fit_Y)
    logORbeta <- V%*%beta
    Q1beta <- as.vector(plogis(qlogis(Q0) + logORbeta ))

  }

  sigma02 <- Q0*(1-Q0)
  sigma1beta2 <- Q1beta*(1-Q1beta)
  omega <- (g0*sigma02 + g1*sigma1beta2)/(g0*sigma02)

  h_star <- V* as.vector(-(g1*sigma1beta2) / (g1*sigma1beta2 + (1-g1)*sigma02))
  H_star <- (A*V + h_star)
  H_star1 <- (V + h_star)
  H_star0 <-  h_star
  offset <- qlogis(Q)
  scale <- apply(V,2, function(v){colMeans_safe(weights* (g1 * sigma1beta2 * v*V))})
  scale_inv <- solve(scale)
  EIFY <- weights*(Delta/G)*omega* (H_star%*%scale_inv) * (Y - Q)



  EIFWA <- apply(V, 2, function(v) {
    (weights*(A*v*(Q1 - Q1beta)) - mean( weights*(A*v*(Q1 - Q1beta))))
  }) %*% scale_inv

  EIF <- EIFY + EIFWA

  var_scaled <- var(EIF)


  fit_Y <- suppressWarnings(glm.fit(A*V,Q, offset = qlogis(Q0), weights = weights,family = binomial(), intercept = F))
  estimates <- coef(fit_Y)
  logORbeta <- V%*%estimates


  compute_predictions <- function(W_newer) {
    V_newer <- model.matrix(formula, data = as.data.frame(W_newer))
    est_grid <-V_newer%*%estimates
    se_grid <- apply(V_newer,1, function(m) {
      sqrt(sum(m * (var_scaled %*%m)))
    } )
    Zvalue <- abs(sqrt(n) * est_grid/se_grid)
    pvalue <- signif(2*(1-pnorm(Zvalue)),5)
    preds_new <- cbind(W_newer,  est_grid,   se_grid/sqrt(n),  se_grid, est_grid - 1.96*se_grid/sqrt(n), est_grid + 1.96*se_grid/sqrt(n),  Zvalue,  pvalue)
    colnames(preds_new) <- c(colnames(W_newer), "estimate",  "se/sqrt(n)", "se", "lower_CI",  "upper_CI", "Z-value", "p-value" )
    return(preds_new)
  }

  preds_new <- compute_predictions(W_new)

  se <- sqrt(diag(var_scaled))
  Zvalue <- abs(sqrt(n) * estimates/se)
  pvalue <- signif(2*(1-pnorm(Zvalue)),5)
  ci <- cbind(estimates,   se/sqrt(n), se ,estimates - 1.96*se/sqrt(n), estimates + 1.96*se/sqrt(n), Zvalue, pvalue)
  colnames(ci) <- c("coefs",   "se/sqrt(n)", "se","lower_CI", "upper_CI", "Z-value", "p-value")
  output <- list(coefs = ci, var_mat = var_scaled, logOR_at_W_new = preds_new, pred_function = compute_predictions,   learner_fits = list(A = fit_A, Y = fit_Y), converged_flag = converged_flag)
  class(output) <- c("spOR","npOddsRatio")
  return(output)
  #######

}




#   compute_predictions <- function(W_newer) {
#     V_newer <- model.matrix(formula, data = as.data.frame(W_newer))
#     est_grid <-V_newer%*%estimates
#     se_grid <- apply(V_newer,1, function(m) {
#       sqrt(sum(m * (var_scaled %*%m)))
#     } )
#     preds_new <- data.frame(W_newer, estimate = est_grid, se = se_grid, lower_CI = est_grid - 1.96*se_grid/sqrt(n), upper_CI = est_grid + 1.96*se_grid/sqrt(n))
#     return(preds_new)
#   }
#
#   preds_new <- compute_predictions(W_new)
#
#   se <- sqrt(diag(var_scaled))
#   ci <- cbind(estimates, se, estimates - 1.96*se/sqrt(n), estimates + 1.96*se/sqrt(n))
#   colnames(ci) <- c("coefs", "se", "lower_CI", "upper_CI")
#   output <- list(coefs = ci, var_mat = var_scaled, logOR_at_W_new = preds_new, pred_function = compute_predictions,   learner_fits = list(A = fit_A, Y = fit_Y), converged_flag = converged_flag)
#   return(output)
#   #######
#
# }
#
