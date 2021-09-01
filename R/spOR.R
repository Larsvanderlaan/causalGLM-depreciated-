

#' Semiparametric Targeted inference for the conditional odds ratio with outcome missingness.
#'
#' Targeted estimates and inference for the odds ratio in a partially-linear logistic-link semiparametric model with `post-treatment`` informative outcome missingness..
#' This version also allows for the outcome is missing-at-random conditional on A,W.
#' The partially-linear logistic model assumes that `logit(P(Y=1|A,W)) = b* A f(W) + h(W)` where `h(W) = logit(P(Y=1|A=0,W))` is unspecified (nonparametric) and `f(W)` is specified by a parametric model.
#' Thus, only a correct parametric model is assumed for the conditional odds ratio, and all other nuisance functions are unspecified (nonparametric).
#'
#' NOTE: For more robust nonparametrically correct inference with no parametric assumptions, use the function \link[npOR] instead. In the \link[npOR] function, the user-specified parametric model is instead treated as an approximation rather than the truth.
#'
#' @param formula An R formula object describing the functional form of the conditional log odds ratio as a fnction of `W`.
#' This corresponds with `f(W)` in the partially linear logistic-link model `logit(P(Y=1|A,W)) = b*Af(W) + h(W)`.
#' @param W A named matrix of baseline covariates
#' @param A A binary vector with values in (0,1) encoding the treatment assignment
#' @param Y A binary outcome variable with values in (0,1)
#' @param Delta A binary vector that takes the value 1 if `Y` is osberved/not-missing and 0 otherwise.
#' @param weights An optional vector of weights for each observation. Use with caution. This can lead to invalid inferences if included naively.
#' @param W_new An optional matrix of new values of the baseline covariates `W` at which to predict odds ratio.
#' @param glm_formula_A (Not recommended). An optional R formula object describing the functional form of P(A=1|W). If provided, \code{glm} is used for the fitting. (Not recommended. This method allows for and works best with flexible machine-learning algorithms.)
#' @param  sl3_learner_A An optional \code{tlverse/sl3} learner object used to estimate P(A=1|W).
#'  If both \code{sl3_learner_A} and \code{glm_formula_A} are not provided, a default learner is used (Lrnr_hal9001).
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
spOR <- function(formula = ~1, W, A, Y, data = NULL, Delta = NULL, weights = NULL, W_new = W, data_new = data, glm_formula_A = NULL, sl3_learner_A = NULL, glm_formula_Y0W = NULL, smoothness_order_Y0W = 1, max_degree_Y0W = 2, num_knots_Y0W = c(25,15,5), reduce_basis = 1e-3, fit_control = list(), sl3_learner_default = Lrnr_hal9001_custom$new(max_degree =2, smoothness_orders = 1, num_knots = c(25,15)), parallel = F,ncores = NULL, targeting_method = c("universal", "iterative"),    fit_Q0W_separate = F, boundsOR = c(1e-3, 1e3), ... ) {
  targeting_method <- match.arg(targeting_method)
  if(parallel) {
    doMC::registerDoMC(ncores)
    fit_control$parallel <- TRUE
  }
  inverse_causal_weighting = FALSE

  W <- as.matrix(W)
  if(is.null(Delta)) {
    Delta <- rep(1, nrow(W))
  }


  formula <- as.formula(formula)
  if(is.null(data)) {
    data <- as.data.frame(W)
  }
  if(is.null(data_new)) {
    data_new <- as.data.frame(W_new)
  }


  V <- model.matrix(formula, data = data)
  V_new <- model.matrix(formula, data = data_new)

  if(is.null(weights)) {
    weights <- rep(1, nrow(W))
  }

  ###### Learn IPCW weights
  # if(inverse_causal_weighting) {
  #   if(sum(1-Delta)>10) {
  #     fit_G <- fit_hal(X = cbind(W,A), Y = Delta, max_degree = 2, family = "binomial", smoothness_orders = 1, num_knots = c(10,2),  fit_control = list(weights=weights))
  #     G <- bound(predict(fit_G, new_data = cbind(W,A)), 0.05)
  #     G <- G[Delta==1]
  #   } else {
  #     G <- 1
  #   }
  #
  # }



  weights <- Delta * weights
  keep <- Delta==1
  if(any(!keep)){
    #Subset to non missing
    V <- V[keep, ,drop = F]
    Delta <- Delta[keep ]
    W_new <- W_new[keep,,drop = F]
    W <- W[keep,,drop = F]
    A <- A[keep ]
    Y <- Y[keep ]
    weights <- weights[keep]

  }
  n <- nrow(W)


  ################################
  #### Learn P(A=1|X) ############
  ################################
  # Default sl3 learner
  if(is.null(sl3_learner_A)){
    sl3_learner_A <- Lrnr_hal9001$new(smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W,   reduce_basis = reduce_basis, ...)
  }
  # If no glm formula then use sl3 learner.
  if(is.null(glm_formula_A)) {

    data_A <- data.frame(W, A)
    covariates_A <- paste0("W", 1:ncol(W))
    colnames(data_A) <- c(covariates_A, "A")
    data_A$weights <- weights
    task_A <- sl3_Task$new(data_A, covariates = covariates_A, outcome = "A", outcome_type = "binomial", weights = "weights" )
    fit_A <- sl3_learner_A$train(task_A)
    g1 <- bound(fit_A$predict(task_A), 0.005)
  } else {
    # use glm is formula supplied
    W_g <- model.matrix(glm_formula_A, data = as.data.frame(W))
    fit_A <- glm.fit(W_g,A, weights = weights, family = binomial(), intercept = F)
    cfs <- coef(fit_A)
    g1 <- bound(as.vector(plogis(W_g%*%cfs)),0.005)
  }




  ################################################################################################
  #### Learn P(Y=1|A,X) under partially linear logistic model assumption #########################
  ################################################################################################
  if(is.null(glm_formula_Y0W)){
    # If no formula use HAL and respect model constraints.

    if(fit_Q0W_separate) {

      fit_control$weights <- weights[A==0]
      fit_Y0 <- fit_hal(X = as.matrix(W)[A==0, ,drop = F],  Y = as.vector(Y)[A==0], family = "binomial", return_x_basis = T, fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W, ...)
      Q0 <- predict(fit_Y0, new_data = as.matrix(W) )
      Q0 <- bound(Q0, 1e-8)
      fit_Y1 <- glm(Y~X-1, data = list(Y=Y[A==1], X=V[A==1, ,drop = F]), weights = weights[A==1], family = binomial(), offset = qlogis(Q0)[A==1])
      beta <- coef(fit_Y1)
      logOR <-  V%*% beta
      Q1 <- plogis(qlogis(Q0) + A*logOR)
      fit_Y <- list(a0 =fit_Y0, a1 =fit_Y1)
    } else {
      if(inverse_causal_weighting) {
        # print("weigting")
        # g <- bound(ifelse(A==1,g1,1-g1),0.05)
        # weight_tmp <- weights*bound(1/g/G,0.05)
      } else {
        weight_tmp <- weights
      }
      fit_control$weights <-weight_tmp
      fit_Y <- fit_hal(X = as.matrix(W), X_unpenalized = as.matrix(A*V), Y = as.vector(Y), family = "binomial", return_x_basis = T, fit_control = fit_control, smoothness_orders = smoothness_order_Y0W, max_degree = max_degree_Y0W, num_knots = num_knots_Y0W, ...)
      Q <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (A*V))
      Q0 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (0*V))
      Q1 <- predict(fit_Y, new_data = as.matrix(W), new_X_unpenalized = (1*V))
      if(inverse_causal_weighting) {
        fit_tmp <- glm.fit(A*V,Y, weights = weights, family = binomial(), offset = qlogis(Q), intercept = F)
        beta <- coef(fit_tmp)
        Q1 <- plogis(qlogis(Q1) + A*V %*% beta)
        Q <- ifelse(A==1,Q1,Q0)
      }

    }


  } else {
    # Use glm if formula supplied
    if(fit_Q0W_separate) {

      W_np <- model.matrix(glm_formula_Y0W, data = as.data.frame(W))
      X <- W_np
      fit_Y0 <- glm.fit(X[A==0,,drop=F],Y[A==0], weights = weights[A==0], family = binomial(), intercept = F)
      cfs <- coef(fit_Y0)
      Q0 <- as.vector(plogis(X%*%cfs))
      fit_Y1 <- glm.fit(V[A==1,,drop=F],Y[A==1], weights = weights[A==1], family = binomial(), intercept = F, offset = qlogis(Q0)[A==1])
      cfs <- coef(fit_Y1)
      Q1 <- as.vector(plogis(qlogis(Q0) + A*V%*%cfs))
      Q <- ifelse(A==1, Q1, Q0)

    } else {

      W_np <- model.matrix(glm_formula_Y0W, data = as.data.frame(W))
      X <- cbind(W_np, A*V)
      X1 <- cbind(W_np, 1*V)
      X0 <- cbind(W_np, 0*V)
      fit_Y <- glm.fit(X,Y, weights = weights, family = binomial(), intercept = F)
      cfs <- coef(fit_Y)
      Q <- as.vector(plogis(X%*%cfs))
      Q1 <- as.vector(plogis(X1%*%cfs))
      Q0 <- as.vector(plogis(X0%*%cfs))

    }
  }


  #### Do some model-compatible bounding of Q0
  denom <- pmax((1-Q1)*(Q0), 1e-8)
  num <- Q1*(1-Q0)
  OR <- num/denom

  logOR <-  log(pmin(pmax(OR, boundsOR[1]), boundsOR[2]) )

  beta <- as.vector(coef(glm(logOR~V-1, family = gaussian(), data = list(V=V, logOR = logOR ))))
  logOR <- V%*%beta
  Q0 <- as.vector(bound(Q0, 1e-4))
  Q1 <- as.vector(plogis(qlogis(Q0) + logOR))
  Q <- ifelse(A==1, Q1, Q0)

  Qinit <- Q
  Qinit1 <- Q1
  Qinit0 <- Q0



  h_star <-  as.vector(-(g1*Qinit1*(1-Qinit1)) / (g1*Qinit1*(1-Qinit1) + (1-g1)*Qinit0*(1-Qinit0)))
  H_star <- V*(A + h_star)
  scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector(Delta * Qinit1*(1-Qinit1) * Qinit0*(1-Qinit0) * g1 * (1-g1) / (g1 * Qinit1*(1-Qinit1) + (1-g1) *Qinit0*(1-Qinit0) )) * v*V)})
  scale_inv <- solve(scale)
  var_unscaled <- as.matrix(var(weights*H_star*(Y-Qinit)))
  var_scaled_init <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)

  ################################
  ##### Targeting Step ###########
  ################################
  converged_flag <- FALSE
  for(i in 1:200) {

    h_star <-  as.vector(-(g1*Q1*(1-Q1)) / (g1*Q1*(1-Q1) + (1-g1)*Q0*(1-Q0)))
    H_star <- V*(A + h_star)
    H_star1 <- V*(1 + h_star)
    H_star0 <- V* h_star
    offset <- qlogis(Q)
    scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector(Delta * Q1*(1-Q1) * Q0*(1-Q0) * g1 * (1-g1) / (g1 * Q1*(1-Q1) + (1-g1) *Q0*(1-Q0) )) * v*V)})

    scale_inv <- solve(scale)

    score <- max(abs(colMeans_safe(weights*H_star%*%scale_inv*as.vector(Y-Q)) ))

    if(abs(score) <= 1/n || abs(score) <= min(0.5,0.5*min(sqrt(diag(var_scaled_init))))/sqrt(n)/log(n)){
      converged_flag <- TRUE
      break
    }
    if(targeting_method == "universal") {
      clev_reduced <- H_star%*%scale_inv
      dir <- colMeans_safe(weights*H_star%*%scale_inv*as.vector(Y-Q))
      dir <- dir/sqrt(mean(dir^2))
      risk <- function(epsilon) {

        Qeps <- plogis(offset + clev_reduced%*%dir * epsilon)
        loss <- -1 * weights * ifelse(Y==1, log(Qeps), log(1-Qeps))
        return(mean(loss))


      }
      optim_fit <- optim(
        par = list(epsilon = 0.01), fn = risk,
        lower = 0, upper = 0.01,
        method = "Brent"
      )
      eps <-  (scale_inv %*%dir) * optim_fit$par

    } else {
      eps <- coef(glm(Y~X-1, family = binomial(), weights =  weights, offset = offset, data = list(Y = Y, X = H_star)))
    }
    #eps <- coef(glm(Y~X-1, family = binomial(), weights = weights, offset = offset, data = list(Y = Y, X = H_star)))

    Q0 <- as.vector(plogis(qlogis(Q0) +  H_star0 %*% eps ))
    Q1 <-  as.vector(plogis(qlogis(Q1) +  H_star1 %*% eps))
    Q <- ifelse(A==1, Q1, Q0)


  }
  if(!converged_flag){

    warning("Targeting did not converge.")
  }



  # Cheap way of extracting coefficients.
  denom <- pmax((1-Q1)*(Q0), 1e-8)
  num <- Q1*(1-Q0)
  OR <- num/denom
  logOR <- log(pmax(OR, 1e-8))
  estimates <- as.vector(coef(glm(logOR~V-1, family = gaussian(), data = list(V=V, logOR = logOR ))))





  h_star <- -(g1*Q1*(1-Q1)) / (g1*Q1*(1-Q1) + (1-g1)*Q0*(1-Q0))
  H_star <- weights*V*(A + h_star)
  scale <- apply(V,2, function(v){colMeans_safe(weights*as.vector(Delta * Q1*(1-Q1) * Q0*(1-Q0) * g1 * (1-g1) / (g1 * Q1*(1-Q1) + (1-g1) *Q0*(1-Q0) )) * v*V)})
  scale_inv <- solve(scale)
  var_unscaled <- as.matrix(var(weights*H_star*(Y-Q)))
  var_scaled <-  scale_inv %*% var_unscaled  %*%  t(scale_inv)


  var_scaled <- var_scaled_init

  compute_predictions <- function(data_new) {
    V_newer <- model.matrix(formula, data = data_new)
    est_grid <-V_newer%*%estimates
    se_grid <- apply(V_newer,1, function(m) {
      sqrt(sum(m * (var_scaled %*%m)))
    } )
    Zvalue <- abs(sqrt(n) * est_grid/se_grid)
    pvalue <- signif(2*(1-pnorm(Zvalue)),5)
    preds_new <- cbind(V_newer,  est_grid,   se_grid/sqrt(n),  se_grid, est_grid - 1.96*se_grid/sqrt(n), est_grid + 1.96*se_grid/sqrt(n),  Zvalue,  pvalue)
    colnames(preds_new) <- c(colnames(V_newer), "estimate",  "se/sqrt(n)", "se", "lower_CI",  "upper_CI", "Z-value", "p-value" )
    return(preds_new)
  }

  preds_new <- compute_predictions(data_new)

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

