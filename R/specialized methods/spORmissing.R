
#' Semiparametric Targeted inference for the conditional odds ratio with post-treatment informed outcome missingness.
#'
#' Semiparametric Targeted estimates and inference for the odds ratio in a partially-linear logistic-link semiparametric model with `post-treatment`` informative outcome missingness.
#' This version also allows for the outcome is missing-at-random conditional on Z,A,W where Z comes after A.
#' The partially-linear logistic model assumes that `logit(P(Y=1|A,W)) = b* A f(W) + h(W)` where `h(W) = logit(P(Y=1|A=0,W))` is unspecified (nonparametric) and `f(W)` is specified by a parametric model.
#' Thus, only a correct parametric model is assumed for the conditional odds ratio, and all other nuisance functions are unspecified (nonparametric).
#'
#' NOTE: For more robust nonparametrically correct inference with no parametric assumptions, use the function \link[npORMissing] instead. In the \link[npORMissing] function, the user-specified parametric model is instead treated as an approximation rather than the truth.
#'
#' @param formula An R formula object describing the functional form of the conditional log odds ratio as a function of `W`.
#' This corresponds with `f(W)` in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`.
#' @param W A named matrix of baseline covariates
#' @param A A binary vector with values in (0,1) encoding the treatment assignment
#' @param Y A binary outcome variable with values in (0,1)
#' @param Z A matrix of post-treatment variables that inform the missingness of the outcome `Y`.
#' This variable can be NULL if the missingness is only effected by `A` and `W`.
#' @param Delta A binary vector that takes the value 1 if `Y` is osberved/not-missing and 0 otherwise.
#' @param weights An optional vector of weights for each observation. Use with caution. This can lead to invalid inferences if included naively.
#' @param W_new An optional matrix of new values of the baseline covariates `W` at which to predict odds ratio.
#' @param glm_formula_A (Not recommended). An optional R formula object describing the functional form of P(A=1|W). If provided, \code{glm} is used for the fitting. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param  sl3_learner_A An optional \code{tlverse/sl3} learner object used to estimate P(A=1|W).
#'  If both \code{sl3_learner_A} and \code{glm_formula_A} are not provided, a default learner is used (Lrnr_hal9001).
#' @param glm_formula_Delta (Not recommended). An optional R formula object describing the functional form of P(Delta=1|A,Z,W) to fit with glm. If provided, it is estimated using glm. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param sl3_learner_Delta An optional \code{tlverse/sl3} learner object used to estimate P(Delta=1|A,Z,W).
#' @param glm_formula_YZ A(Not recommended). An optional R formula object describing the functional form of the full conditional outcome mean P(Y=1|Z,A,W). If provided, it is estimated using glm. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param sl3_learner_YZ An optional \code{tlverse/sl3} learner object used to estimate the full conditional outcome mean P(Y=1|Z,A,W).
#' @param glm_formula_Y0W (Not recommended). An optional R formula object describing the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`
#' @param smoothness_order_Y0W Smoothness order of the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)` to be estimated nonparametrically using the Highly Adaptive Lasso (\link{hal9001}), a powerful spline regression algorithm. 0 = discontinuous piece-wise constant function, 1 = continuous piece-wise linear, 2 = smooth piece-wise quadratic
#' @param max_degree_Y0W Max degree of interaction (of spline basis functions) of the nuisance function `h(W) := logit(P(Y=1|A=0,W))`in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)` to be estimated nonparametrically using the Highly Adaptive Lasso (\link{hal9001}). `max_degree=1` corresponds with an additive model, `max_degree=2` corresponds with a bi-additive (two-way) model. This parameter significantly affects computation time.
#' @param num_knots_Y0W A vector specifying the number of knots to use when generating `HAL` spline basis functions of each interaction degree. For computational benefits, the number of knots should decrease exponentially with degree.
#' @param reduce_basis See analagous argument in package \link{hal9001}.
#' @param fit_control See analagous argument in package \link{hal9001}.
#' @param ... Other arguments to be passed to \link{hal9001::fit_hal} for fitting.
#' @param sl3_learner_default A default sl3 Learner to be used if neither a glm formula or sl3 learner is provided for one of the nuisance functions.
#' By default, Lrnr_hal9001 is used.
#' @importFrom doMC registerDoMC
#' @export
spORmissing <- function(formula = logOR~1, W, A, Y, Z= stop("No Z variable given. Use spOR or npOR instead if the variable `Z` is not used"), Delta, weights = NULL, W_new = W, glm_formula_A = NULL, sl3_learner_A = NULL, glm_formula_Delta= NULL, sl3_learner_Delta = NULL,  glm_formula_YZ = NULL, sl3_learner_YZ = NULL, glm_formula_Y0W = NULL, smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(20,5), reduce_basis = 1e-3, fit_control = list(), sl3_learner_default = Lrnr_hal9001_custom$new(max_degree =2, smoothness_orders = 1, num_knots = c(30,10)), parallel = F,ncores = NULL, ... ) {

  if(parallel) {
    doMC::registerDoMC(ncores)
    fit_control$parallel <- TRUE
  }

  glm_formula_Y_W <- glm_formula_Y0W
  W <- as.matrix(W)
  Z <- as.matrix(Z)
  if(is.null(Delta)) {
    Delta <- rep(1, nrow(W))
  }

  n <- nrow(W)
  formula <- as.formula(formula)
  V <- model.matrix(formula, data = as.data.frame(W))
  V_new <- model.matrix(formula, data = as.data.frame(W_new))
  if(is.null(weights)) {
    weights <- rep(1, nrow(W))
  }

  keep <- Delta==1

  ####################################
  ############ Default learners ############
  ####################################
  if(is.null(sl3_learner_A)){
    sl3_learner_A <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W,   reduce_basis = reduce_basis, ...)
  }
  if(is.null(sl3_learner_Delta)){
    sl3_learner_Delta <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W,   reduce_basis = reduce_basis, ...)
  }
  if(is.null(sl3_learner_YZ)){
    sl3_learner_YZ <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W,   reduce_basis = reduce_basis, ...)
  }

  ################################################################################################
  #### Learn P(Y=1|Z, A,W)  #########################
  ################################################################################################
  print("Learning barQ")
  if(is.null(glm_formula_YZ)){
    # If no formula use HAL and respect model constraints.
    data_YZ <- data.frame(W, A, Z)
    covariates_YZ <- c(paste0("W", 1:ncol(W)), "A",paste0("Z", 1:ncol(Z)))
    colnames(data_YZ) <- c(covariates_YZ)
    data_YZ$Y <- Y
    data_YZ$weights <- weights
    task_YZ <- sl3_Task$new(data_YZ, covariates = covariates_YZ, outcome = "Y", outcome_type = "binomial", weights = "weights" )
    data_YZ0 <- data_YZ
    data_YZ0$A <- 0
    task_YZ0 <- sl3_Task$new(data_YZ0, covariates = covariates_YZ, outcome = "Y", outcome_type = "binomial", weights = "weights" )
    data_YZ1 <- data_YZ
    data_YZ1$A <- 1
    task_YZ1 <- sl3_Task$new(data_YZ1, covariates = covariates_YZ, outcome = "Y", outcome_type = "binomial", weights = "weights" )
    fit_YZ <- sl3_learner_YZ$train(task_YZ[keep])
    barQ <- bound(fit_YZ$predict(task_YZ),0.0005)
    print(quantile(barQ))
    barQ1 <- bound(fit_YZ$predict(task_YZ1),0.0005)
    barQ0 <- bound(fit_YZ$predict(task_YZ0),0.0005)
  } else {
    # Use glm if formula supplied

    W_np <- model.matrix(glm_formula_YZ, data = as.data.frame(cbind(W,Z)))
    X <- cbind(W_np, A*W_np)
    X1 <- cbind(W_np, 1*W_np)
    X0 <- cbind(W_np, 0*W_np)
    fit_Y <- glm.fit(X,Y, weights = Delta*weights, family = binomial(), intercept = F)
    cfs <- coef(fit_Y)
    barQ <- bound(as.vector(plogis(X%*%cfs)),0.0005)
    barQ1 <- bound(as.vector(plogis(X1%*%cfs)),0.0005)
    barQ0 <- bound(as.vector(plogis(X0%*%cfs)),0.0005)
  }

  print("Learning Q")
  ################################################################################################
  #### Learn P(Y=1|A,X) under partially linear logistic model assumption #########################
  ################################################################################################
  if(is.null(glm_formula_Y_W)){
    # If no formula use HAL and respect model constraints.
    fit_control$weights <- weights
    fit_Y <- suppressWarnings(fit_hal(X = as.matrix(W), X_unpenalized = as.matrix(A*V),
                                      Y = as.vector(barQ), family = binomial(), fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W, ...))
    Q <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = as.matrix(A*V))
    Q0 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = as.matrix(0*V))
    Q1 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = as.matrix(1*V))

  } else {
    # Use glm if formula supplied
    W_np <- model.matrix(glm_formula_Y_W, data = as.data.frame(W))
    X <- cbind(W_np, A*V)
    X1 <- cbind(W_np, 1*V)
    X0 <- cbind(W_np, 0*V)
    fit_Y <- glm.fit(X,barQ, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_Y)
    Q <- as.vector(plogis(X%*%cfs))
    Q1 <- as.vector(plogis(X1%*%cfs))
    Q0 <- as.vector(plogis(X0%*%cfs))
  }
  denom <- pmax((1-Q1)*(Q0), 1e-8)
  num <- Q1*(1-Q0)
  OR <- num/denom
  logOR <-  log(pmin(pmax(OR, 1e-8), 5000) )
  beta <- as.vector(coef(glm(logOR~V-1, family = gaussian(), data = list(V=V, logOR = logOR ))))
  logOR <- V%*%beta
  Q0 <- bound(Q0, 1e-6)
  Q1 <- plogis(qlogis(Q0) + logOR)
  Q <- ifelse(A==1, Q1, Q0)
  ################################
  #### Learn P(A=1|X) ############
  ################################

  # If no glm formula then use sl3 learner.
  if(is.null(glm_formula_Delta)) {

    data_Delta <- data.frame(W, A, Z)
    covariates_Delta  <- c(paste0("W", 1:ncol(W)), "A",paste0("Z", 1:ncol(Z)))
    colnames(data_Delta ) <- covariates_Delta
    data_Delta$Delta <- Delta
    data_Delta$weights <- weights
    task_Delta <- sl3_Task$new(data_Delta, covariates = covariates_Delta, outcome = "Delta", outcome_type = "binomial", weights = "weights" )
    fit_Delta <- sl3_learner_Delta$train(task_Delta)
    G <- fit_Delta$predict(task_Delta)
    print(range(G))
    G <- pmax(G, 0.005)
  } else {
    # use glm is formula supplied
    W_np <- model.matrix(glm_formula_Delta, data = as.data.frame(cbind(W,Z)))
    X <- cbind(W_np, A*W_np)
    X1 <- cbind(W_np, 1*W_np)
    X0 <- cbind(W_np, 0*W_np)
    fit_Delta<- glm.fit(X,Delta, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_Delta)
    G <- as.vector(plogis(X%*%cfs))
    G <- pmax(G, 0.005)
  }
  if(is.null(glm_formula_A)) {

    data_A <- data.frame(W, A)
    covariates_A <- paste0("W", 1:ncol(W))
    colnames(data_A) <- c(covariates_A, "A")
    data_A$weights <- weights
    task_A <- sl3_Task$new(data_A, covariates = covariates_A, outcome = "A", outcome_type = "binomial", weights = "weights" )
    fit_A <- sl3_learner_A$train(task_A)
    g1 <- bound(fit_A$predict(task_A),0.0005)
  } else {
    # use glm is formula supplied
    W_g <- model.matrix(glm_formula_A, data = as.data.frame(W))
    fit_A <- glm.fit(W_g,A, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_A)
    g1 <- bound(as.vector(plogis(W_g%*%cfs)),0.0005)
  }






  ################################
  ##### Targeting Step ###########
  ################################
  converged_flag <- FALSE
  for(i in 1:50) {


    ################################
    ##### Update  Q ###########
    ################################
    for(j in 1:10) {

      h_star <-  as.vector(-(g1*Q1*(1-Q1)) / (g1*Q1*(1-Q1) + (1-g1)*Q0*(1-Q0)))
      H_star <- V*(A + h_star)
      H_star1 <- V*(1 + h_star)
      H_star0 <-  V*h_star
      scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector( Q1*(1-Q1) * Q0*(1-Q0) * g1 * (1-g1) / (g1 * Q1*(1-Q1) + (1-g1) *Q0*(1-Q0) )) * v*V)})
      scale_inv <- solve(scale)
      var_unscaled <- as.matrix(var(weights*H_star*(Y-Q)))
      var_scaled <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)

      score <- sum(abs(colMeans_safe(weights*H_star%*%scale_inv*as.vector(barQ-Q)) ))
      print("Inner")
      print(score)
      if(abs(score) <= min(0.5,0.5*mean(sqrt(diag(var_scaled))))/sqrt(n)/log(n)){
        break
      }
      offset <- qlogis(Q)
      eps <- coef(suppressWarnings(glm(Y~X-1, family = binomial(), weights = weights, offset = offset, data = list(Y = barQ, X = H_star))))
      Q <- as.vector(plogis(offset +  H_star %*% eps))
      Q0 <- as.vector(plogis(qlogis(Q0) +  H_star0 %*% eps ))
      Q1 <- as.vector(plogis(qlogis(Q1) +  H_star1 %*% eps))
    }

    ################################
    ##### Update barQ ###########
    ################################
    h_star <- as.vector(-(g1*Q1*(1-Q1)) / (g1*Q1*(1-Q1) + (1-g1)*Q0*(1-Q0)))
    H_star <- V*(A + h_star)
    H_star1 <- V*(1 + h_star)
    H_star0 <-  V*h_star

    scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector( Q1*(1-Q1) * Q0*(1-Q0) * g1 * (1-g1) / (g1 * Q1*(1-Q1) + (1-g1) *Q0*(1-Q0) )) * v*V)})
    scale_inv <- solve(scale)
    var_unscaled <- as.matrix(var(weights*H_star*(Y-Q)))
    var_scaled <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)

    offset <- qlogis(barQ)
    eps <- coef(suppressWarnings(glm(Y~X-1, family = binomial(), weights = Delta*weights/G, offset = offset, data = list(Y = Y, X = H_star))))
    barQ <- as.vector(plogis(offset +  H_star %*% eps))
    barQ0 <- as.vector(plogis(qlogis(barQ0) +  H_star0 %*% eps ))
    barQ1 <- as.vector(plogis(qlogis(barQ1) +  H_star1 %*% eps))


    ################################
    ##### Check convergence ###########
    ################################
    h_star <- -1* as.vector((g1*Q1*(1-Q1)) / (g1*Q1*(1-Q1) + (1-g1)*Q0*(1-Q0)))
    H_star <- V*(A + h_star)
    H_star1 <- V*(1 + h_star)
    H_star0 <-  V*h_star

    scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector(Q1*(1-Q1) * Q0*(1-Q0) * g1 * (1-g1) / (g1 * Q1*(1-Q1) + (1-g1) *Q0*(1-Q0) )) * v*V)})
    scale_inv <- solve(scale)
    EIF1 <- (Delta/G)*weights*H_star%*%scale_inv*as.vector(Y-barQ)
    EIF2 <- weights*H_star%*%scale_inv*as.vector(barQ-Q)
    EIF <- EIF1 + EIF2

    var_scaled <- as.matrix(var(EIF))
    #var_scaled <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)

    score <- sum(abs(colMeans_safe(EIF1 + EIF2 )))

    print(score)
    if(score <= min(0.5,0.5*mean(sqrt(diag(var_scaled))))/sqrt(n)/log(n)) {
      print("converged")
      converged_flag <- TRUE
      break
    }



  }


  # Cheap way of extracting coefficients.
  denom <- pmax((1-Q1)*(Q0), 1e-8)
  num <- Q1*(1-Q0)
  OR <- num/denom
  logOR <- log(pmax(OR, 1e-8))
  estimates <- as.vector(coef(suppressWarnings(glm(logOR~V-1, family = gaussian(), data = list(V=V, logOR = logOR )))))



  var_scaled <-  var_scaled

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
  class(output) <- c("spORMissing","npOddsRatio")
  return(output)
  #######

}



  # compute_predictions <- function(W_newer) {
  #   V_newer <- model.matrix(formula, data = as.data.frame(W_newer))
  #   est_grid <-V_newer%*%estimates
  #   se_grid <- apply(V_newer,1, function(m) {
  #     sqrt(sum(m * (var_scaled %*%m)))
  #   } )
  #   preds_new <- data.frame(W_newer, estimate = est_grid, se = se_grid, lower_CI = est_grid - 1.96*se_grid/sqrt(n), upper_CI = est_grid + 1.96*se_grid/sqrt(n))
  #   return(preds_new)
  # }
  #
  # preds_new <- compute_predictions(W_new)
  #
  # se <- sqrt(diag(var_scaled))
  # ci <- cbind(estimates, se, estimates - 1.96*se/sqrt(n), estimates + 1.96*se/sqrt(n))
  # colnames(ci) <- c("coefs", "se", "lower_CI", "upper_CI")
  # output <- list(coefs = ci, var_mat = var_scaled, logOR_at_W_new = preds_new, pred_function = compute_predictions,   learner_fits = list(A = fit_A, Y = fit_Y), converged_flag = converged_flag)
  # return(output)
  #######


