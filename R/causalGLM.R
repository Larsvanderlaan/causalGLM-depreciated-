
#' causalGLM
#' A default and more user-friendly front-end of the implemented methods.
#' Supports conditional average treatment effect (CATE), conditional odds ratio (OR), and conditional relative risk (RR) estimation
#' Highly Adaptive Lasso (HAL, R package: tlverse/hal9001), a flexible and adaptive spline regression estimator, is used as default learner (piece-wise linear).
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param W A named matrix or data.frame of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y n outcome variable (continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "autoHAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001} Good for most sample sizes when propertly tuned. See arguments \code{max_degree_Y0W} and \code{num_knots_Y0W}.
#' "glm": Fit nuisances with parametric model. Best for smaller sample sizes (e.g. n =30-100). See arguments \code{glm_formula_A}, \code{glm_formula_Y} and \code{glm_formula_Y0}.
#' "glmnet": Learn using lasso with glmnet. Best for smaller sample sizes (e.g. n =30-100)
#' "gam": Learn using generalized additive models with mgcv. Good for small-to-medium-small sample sizes.
#' "mars": Multivariate adaptive regression splines with \code{earth}. Good for small-to-medium-small sample sizes.
#' "ranger": Robust random-forests with the package \code{Ranger} Good for medium-to-large sample sizes.
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7. Good for medium-to-large sample sizes.
#' We recommend performing simulations checking 95% CI coverage when choosing learners (especially in smaller sample sizes).
#' @param pool_A_when_training Default: TRUE. This argument is ignored if \code{learning_method} = `autoHAL`, `glm`, `gam`, or `glmnet` and \code{sl3_Learner_Y0W} is NULL. 
#' Otherwise this is a boolean for whether to estimate the conditional mean/regression of Y by combining observations with A=0,A=1 ...
#' Or to estimate E[Y|A=1,W] and E[Y|A=0,W] with separate regressions (this is nonparametric in the interaction with A).
#' When \code{pool_A_when_training} is TRUE, the design matrix passed to the regression algorithm/learner is cbind(W,A*V) where V is the design matrix specified by the argument \code{formula}.
#' Therefore, it may not be necessary to use learners that model (treatment) interactions when this argument is TRUE.
#' For \code{learning_method} = glm, gam,  mars, and glmnet this argument is set to TRUE automatically.
#' In high dimensions, pool_A_when_training = FALSE may be preferred to prevent dilution of the treatment interactions in the fitting.
#' @param estimand Estimand/parameter to estimate. Choices are:
#' CATE: Estimate conditional average treatment effect with \code{spCATE} assuming it satisfies parametric model \code{formula}.
#' OR: Estimate conditional odds ratio with \code{spOR} assuming it satisfies parametric model \code{formula}.
#' OR: Estimate conditional relative risk with \code{spRR} assuming it satisfies parametric model \code{formula}.
#' @param cross_fit Whether to cross-fit the initial estimator. This is always set to FALSE if argument \code{sl3_Learner} is provided.
#'  learning_method = AutoML (default) does use cross-fitting for computational reasons and due to it being unnecessary for HAL.
#'  learning_method = `xgboost` and `ranger` are always cross-fitted regardless of the value of \code{cross_fit}
#'  All other learning_methods are only cross-fitted if `cross_fit=TRUE`. 
#'  Note, it is not necessary to cross-fit glm, glmnet, gam or mars as long as the dimension of W is not too high.
#'  In smaller samples and lower dimensions, it may fact hurt to cross-fit.  
#' @param sl3_Learner_A A \code{sl3} Learner object to use to estimate nuisance function P(A=1|W) with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided. 
#' If you wish to cross-fit the learner \code{sl3_Learner} then do: sl3_Learner <- Lrnr_cv$new(sl3_Learner).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param sl3_Learner_Y A \code{sl3} Learner object to use to estimate nuisance functions [Y|A=1,W] and E[Y|A=0,W] (depending on method) with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided. 
#' Keep in mind the value of the argument \code{pool_A_when_training}. If FALSE  then E[Y|A=0,W] is estimated by itself. 
#' Therefore, it may not be needed to add interactions, since treatment interactions are automatic by stratification.
#'If TRUE, the design matrix passed to the pooled learner contains A*V where V is the design matrix obtained from \code{formula}.
#' For some learners, it may also be unnecessary to include interactions in this case.
#' #' If you wish to cross-fit the learner \code{sl3_Learner} then do: sl3_Learner <- Lrnr_cv$new(sl3_Learner).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param glm_formula_A A glm formula for P(A=1|W). Only used if learning_method = "glm".
#' @param glm_formula_Y A glm formula for E[Y|A,W] or E[Y|A=0,W] depending on method. Only used if learning_method = "glm".
#' @param weights An optional vector of weights to use in procedure.
#' @param parallel Whether to parallize HAL whenever it is used in fitting.
#' @param ncores Number of cores to use in parallelization.
#' @param smoothness_order_Y0W Smoothness order for HAL estimator of E[Y|A=0,W] (see \code{hal9001/fit_hal})
#' smoothness_order_Y0W = 1 is piece-wise linear. smoothness_order_Y0W = 0 is piece-wise constant.
#' @param max_degree_Y0W Max interaction degree for HAL estimator of E[Y|A=0,W] (see \code{hal9001/fit_hal})
#' @param num_knots_Y0W Number of knots by interaction degree for HAL estimator of E[Y|A=0,W](see \code{hal9001/fit_hal}). Used to generate basis functions.
#' num_knots_Y0W = c(1) is equivalent to main term glmnet/LASSO. (Assuming max_degree_Y0W = 1)
#' num_knots_Y0W = c(1,1) is equivalent to glmnet/LASSO with both main-terms and all two-way interactions (e.g. Y~ W1 + W1 + W1*W2 + ...).  (Assuming max_degree_Y0W = 2)
#' num_knots_Y0W = c(10) is an additive piece-wise linear model with 10 knot points.  (Assuming max_degree_Y0W = 1)
#' num_knots_Y0W = c(10,5) is a bi-additive model with the same one-way basis functions as above, but also two-way interaction piece-wise linear basis functions generated by main-term one-way basis functions with 5 knots. (Assuming max_degree_Y0W = 2)
#' num_knots_Y0W = c(10,5,1) generates same basis functions as above and also the three-way interaction basis functions with only a single knot point at the origin (e.g. the triple interaction `W1*W2*W3`) (Assuming max_degree_Y0W = 3)
#' 
#' @param data_list A named list containing the arguments `W`, `A` and `Y`. For example, data_list = list(W = data[,c("W1", "W2")], A = data[,"A"], Y = data[,"Y"])
#' @param ... Other arguments to pass to main routine (spCATE, spOR, spRR) 
#' @export
causalGLM <- function(formula, W, A, Y, estimand = c("CATE", "OR", "RR"),   learning_method = c("autoHAL", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"), pool_A_when_training = TRUE,  cross_fit = ifelse(ncol(W) >= 12, T, F),  sl3_Learner_A = NULL, sl3_Learner_Y = NULL, glm_formula_A = NULL, glm_formula_Y = NULL,  weights = NULL, data_list = NULL, parallel =  F, ncores = NULL, smoothness_order_Y0W = 1, max_degree_Y0W = ifelse(nrow(W) >= 200, 2,1), num_knots_Y0W = c(ifelse(nrow(W) >= 500 && ncol(W) <= 20, 20, 10),5,1), constant_variance_CATE = FALSE, ... ){
  smoothness_order_Y0W <- smoothness_order_Y0W[1]
  fast_analysis = TRUE
  if(!(smoothness_order_Y0W %in% c(0,1))) {
    stop("smoothness_order_Y0W must be 0 or 1.")
  }
  estimand <- match.arg(estimand)
  learning_method <- match.arg(learning_method)
  family <- NULL
  if(estimand == "RR") {
    family <- poisson()
    if(learning_method %in% c("glmnet", "autoHAL")) {
      family <- "poisson"
    }
  }
  
  num_knots_Y0W <- num_knots_Y0W[1:max_degree_Y0W]
   
  if(learning_method %in% c("glm", "glmnet", "gam","mars")) {
    pool_A_when_training <- TRUE
  }
  inference_type =  "semiparametric"
  if(parallel & !is.null(ncores)) {
    doMC::registerDoMC(ncores)
  } 
  if(!is.null(data_list)) {
    W <- data_list$W
    A <- data_list$A 
    Y <- data_list$Y
  }
  
  if( inference_type=="parametric") {
    stop("`parametric`` is currently not supported.  You can set learning_method = `glm` to get nearly the same estimates as parametric glm, although this is not recommended (use MARS or glmnet instead for adaptive parametric fitting).")
    #warning("If inference_type is `parametric` then argument `formula` should be a full formula with both A and Y. I will just return glm output.")
  } 
  W <- as.matrix(W)
  n <- nrow(W)
  p <- ncol(W)
  sl3_Learner <- NULL
  sl3_Learner_Y_orig <- sl3_Learner_Y
  sl3_Learner_A_orig <- sl3_Learner_A
 
  if(is.null(sl3_Learner_Y) || is.null(sl3_Learner_A)) {
    if(learning_method == "autoHAL" ) {
      sl3_Learner_A <- autoML(n,p, parallel, fast_analysis, NULL)
      sl3_Learner_Y <- NULL #autoML(n,p, parallel, fast_analysis, family)
    } else if(learning_method == "glmnet" ) {
      sl3_Learner_A <- Lrnr_glmnet$new()
      sl3_Learner_Y <- Lrnr_glmnet$new(family = family)
      
    } else if(learning_method == "glm" ) {
      sl3_Learner_A <- Lrnr_glm$new(formula = glm_formula_A)
      sl3_Learner_Y <- Lrnr_glm$new(formula = glm_formula_Y, family = family)
      
      
    } else if(learning_method == "gam" ) {
      sl3_Learner_A <- Lrnr_gam$new()
      sl3_Learner_Y <- Lrnr_gam$new(family = family)
    } else if(learning_method == "mars" ) {
      sl3_Learner <- Lrnr_earth$new()
    } else if (learning_method == "ranger") {
      sl3_Learner <- Lrnr_cv$new(Lrnr_ranger$new())
    } else if(learning_method == "xgboost" ) {
      objective <- NULL
      if(estimand == "RR") {
        objective <- "count:poisson"
      }
      sl3_Learner <- Stack$new( Lrnr_glmnet$new(family = family), Lrnr_xgboost$new(max_depth =3,  objective = objective), Lrnr_xgboost$new(max_depth =4, objective = objective), Lrnr_xgboost$new(max_depth =5, objective = objective), Lrnr_xgboost$new(max_depth =6, objective = objective), Lrnr_xgboost$new(max_depth =7, objective = objective))
      sl3_Learner_Y <- make_learner(Pipeline, Lrnr_cv$new(sl3_Learner), Lrnr_cv_selector$new(loss_squared_error))
      sl3_Learner <- Stack$new( Lrnr_glmnet$new(family = NULL), Lrnr_xgboost$new(max_depth =3,  objective = NULL), Lrnr_xgboost$new(max_depth =4, objective = NULL), Lrnr_xgboost$new(max_depth =5, objective = NULL), Lrnr_xgboost$new(max_depth =6, objective = NULL), Lrnr_xgboost$new(max_depth =7, objective = NULL))
      sl3_Learner_A <- make_learner(Pipeline, Lrnr_cv$new(sl3_Learner), Lrnr_cv_selector$new(loss_squared_error))
      
    }
    if(is.null(sl3_Learner_Y)) {
      sl3_Learner_Y <- sl3_Learner
    }
    if( is.null(sl3_Learner_A)) {
      sl3_Learner_A <- sl3_Learner
    }
    if(cross_fit & learning_method %in% c("glm", "gam", "glmnet", "mars") ) {
      sl3_Learner_A <- Lrnr_cv$new(sl3_Learner_A)
      sl3_Learner_Y <- Lrnr_cv$new(sl3_Learner_Y)
    }
    if(!is.null(sl3_Learner_Y_orig)) {
      sl3_Learner_Y <- sl3_Learner_Y_orig
    }
    if( !is.null(sl3_Learner_A_orig)) {
      sl3_Learner_A <- sl3_Learner_A_orig
    }
    
    
  }
  
  sl3_Learner_Y0W <- sl3_Learner_Y
 
  if(estimand == "RR") {
    return(spRR(formula_logRR =  formula, W, A, Y, pool_A_when_training = pool_A_when_training, sl3_Learner_A = sl3_Learner_A, sl3_Learner_Y = sl3_Learner_Y,   weights = weights,  smoothness_order_Y0W = smoothness_order_Y0W, max_degree_Y0W = max_degree_Y0W, num_knots_Y0W = num_knots_Y0W,  fit_control = list(parallel = parallel),...))
  }
  if(inference_type == "semiparametric") {
    if(estimand == "CATE") {
      return(spCATE(formula_CATE =  formula, W, A, Y, pool_A_when_training = pool_A_when_training,  sl3_Learner_A = sl3_Learner_A, sl3_Learner_Y = sl3_Learner_Y,   weights = weights,  smoothness_order_Y0W = smoothness_order_Y0W, max_degree_Y0W = max_degree_Y0W, num_knots_Y0W = num_knots_Y0W,  fit_control = list(parallel = parallel), constant_variance = constant_variance_CATE , ...))
    }
    if(estimand == "OR") {
      return(spOR(formula_logOR = formula, W, A, Y,  pool_A_when_training = pool_A_when_training, weights = weights, W_new = W,  sl3_Learner_A = sl3_Learner_A, sl3_Learner_Y0W = sl3_Learner_Y0W,  glm_formula_Y0W = glm_formula_Y0W, smoothness_order_Y0W = smoothness_order_Y0W, max_degree_Y0W = max_degree_Y0W, num_knots_Y0W = num_knots_Y0W, reduce_basis = 1e-3, fit_control = list(parallel = parallel), sl3_learner_default = sl3_Learner_Y0W ,... ))
    }
  } 
  
  
  
}




#' causalGLMwithLASSO
#' causalGLM in high dimensions. Valid inference with data-adaptive variable selection.
#' A wrapper for causalGLM with partially-penalized \code{glmnet}/LASSO as base learner (Interactions terms with `A` are not penalized). 
#' This method is useful for high dimensional settings where other learners are slow or poorly behaved.
#' It may also be useful for smaller sample sizes where other machine-learning algorithms may overfit.
#' Otherwise, we do not recommend using this function for lower dimensional settings since glmnet can be mispecified.
#' 
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param W A named matrix or data.frame of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y n outcome variable (continuous, nonnegative or binary depending on method)
#' @param estimand Estimand/parameter to estimate. Choices are:
#' CATE: Estimate conditional average treatment effect with \code{spCATE} assuming it satisfies parametric model \code{formula}.
#' OR: Estimate conditional odds ratio with \code{spOR} assuming it satisfies parametric model \code{formula}.
#' OR: Estimate conditional relative risk with \code{spRR} assuming it satisfies parametric model \code{formula}.
#' @param cross_fit Whether to cross-fit the initial estimator. By default, TRUE. 
#' In lower dimensions, we recommend setting this to FALSE.
#' @param weights An optional vector of weights to use in procedure.
#' @param data_list A named list containing the arguments `W`, `A` and `Y`. For example, data_list = list(W = data[,c("W1", "W2")], A = data[,"A"], Y = data[,"Y"])
#' @param ... Other arguments to pass to glmnet (NOTE: this use is different than that of \code{causalGLM})
#' @export
causalGLMwithLASSO <- function(formula, W, A, Y, estimand = c("CATE", "OR", "RR"), cross_fit = TRUE,weights = NULL,data_list = NULL, constant_variance_CATE = FALSE,  return_competitor  = F,...  )  {
  V <- model.matrix(formula, as.data.frame(W))
  penalty.factor <- c(rep(1, ncol(W)), rep(1e-10, ncol(V)))
  
  if(estimand == "RR") {
    family <- "poisson"
  } else {
    family <- NULL
  }
  lrnr <- Lrnr_glmnet$new(family = family, penalty.factor=penalty.factor,...)
  if(cross_fit) {
    lrnr <- Lrnr_cv$new(lrnr)
  }
  causalGLM(formula, W, A, Y, estimand,sl3_Learner_Y = lrnr, learning_method = "glmnet", cross_fit = cross_fit, weights = weights, num_knots_Y0W = 1, max_degree_Y0W =1, data_list = data_list , constant_variance = constant_variance_CATE, return_competitor = return_competitor )
}



 

#'  
autoML <- function(n, p, parallel = F, fast_analysis = F, family = NULL) {
  fit_control <- list(parallel = parallel)
  
  if(p >= 200) {
    lrnr <- Lrnr_glmnet$new()
    return(lrnr)
  }
  
  if(p >= 50) {
    lrnr <-  Lrnr_hal9001$new(family=family, smoothness_orders = 1, max_degree = 1, num_knots =  c(5,0), fit_control = fit_control)
    return(lrnr)
  }
  
  if(p >= 25) {
    max_degree <- 1
    
  } else {
    max_degree <- 2
    
  }
  
  if(n<=50) {
    lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = 1, num_knots =  c(1,0), fit_control = fit_control)
  }
  else if(n<=100) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = 2, num_knots =  c(1,1), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = 1, num_knots =  c(1,1), fit_control = fit_control)
    }
  } else if(n<=250) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = 2, num_knots =  c(7,3), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = max_degree, num_knots =  c(3,1), fit_control = fit_control)
    }
  } else if(n<=500) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = 2, num_knots =  c(12,5), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = max_degree, num_knots =  c(10,5), fit_control = fit_control)
    }
  } else if(n<=1000) {
    lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = max_degree, num_knots =  c(20,10), fit_control = fit_control)
  } else {
    if(fast_analysis) {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = max_degree, num_knots =  c(10,5), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(family=family,smoothness_orders = 1, max_degree = max_degree, num_knots =  c(25,15), fit_control = fit_control)
    }
  }
  return(lrnr)
}



