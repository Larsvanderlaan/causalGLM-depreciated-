 
#' causalGLM
#' A default and more user-friendly front-end of the implemented methods.
#' 
#' @param formula A R formula object specifying the parametric form of CATE, OR, or RR (depending on method).
#' @param W A named matrix or data.frame of baseline covariates to condition on.
#' @param A A binary treatment assignment vector
#' @param Y n outcome variable (continuous, nonnegative or binary depending on method)
#' @param learning_method Machine-learning method to use. This is overrided if argument \code{sl3_Learner} is provided. Options are:
#' "autoHAL": Adaptive robust automatic machine-learning using the Highly Adaptive Lasso \code{hal9001}
#' "glm": Fit nuisances with parametric model. See arguments \code{glm_formula_A}, \code{glm_formula_Y} and \code{glm_formula_Y0}.
#' "glmnet": Learn using lasso with glmnet.
#' "gam": Learn using generalized additive models with mgcv.
#' "mars": Multivariate adaptive regression splines with \code{earth}.
#' "ranger": Robust random-forests with the package \code{Ranger}
#' "xgboost": Learn using a default cross-validation tuned xgboost library with max_depths 3 to 7.
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
#' @param sl3_Learner A \code{sl3} Learner object to use to estimate nuisance functions with machine-learning.
#' Note, \code{cross_fit} is automatically set to FALSE if this argument is provided. 
#' If you wish to cross-fit the learner \code{sl3_Learner} then do: sl3_Learner <- Lrnr_cv$new(sl3_Learner).
#' Cross-fitting is recommended for all tree-based algorithms like random-forests and gradient-boosting.
#' @param glm_formula_A A glm formula for P(A=1|W). Only used if learning_method = "glm".
#' @param glm_formula_Y A glm formula for E[Y|A,W]. Only used if learning_method = "glm".
#' @param glm_formula_Y0W A glm formula for E[Y|A=0,W]. Only used if learning_method = "glm".
#' @param weights An optional vector of weights to use in procedure.
#' @param parallel Whether to parallize HAL whenever it is used in fitting.
#' @param ncores Number of cores to use in parallelization.
#' @param smoothness_order Smoothness order for HAL (see \code{hal9001/fit_hal})
#' @param max_degree Max interaction degree for HAL (see \code{hal9001/fit_hal})
#' @param num_knots Number of knots by interaction degree for HAL (see \code{hal9001/fit_hal}). Used to generate basis functions.
#' @export
causalGLM <- function(formula, W, A, Y, estimand = c("CATE", "OR", "RR"),   learning_method = c("autoHAL", "glm", "glmnet", "gam", "mars", "ranger", "xgboost"),    cross_fit = ifelse(ncol(W) >= 12, T, F),  sl3_Learner = NULL, glm_formula_A = NULL, glm_formula_Y = NULL, glm_formula_Y0W = glm_formula_Y, weights = NULL, data_list = NULL,  fast_analysis = TRUE, parallel =  F, ncores = NULL, smoothness_order = 1, max_degree = 2, num_knots = c(10,5) ){
  inference_type =  "semiparametric"
  if(parallel & !is.null(ncores)) {
    doMC::registerDoMC(ncores)
  } 
  if(!is.null(data_list)) {
    W <- data_list$W
    A <- data_list$A 
    Y <- data_list$Y
  }
  estimand <- match.arg(estimand)
  
  learning_method <- match.arg(learning_method)
  if( inference_type=="parametric") {
    stop("`parametric`` is currently not supported.  You can set learning_method = `glm` to get nearly the same estimates as parametric glm, although this is not recommended (use MARS or glmnet instead for adaptive parametric fitting).")
    #warning("If inference_type is `parametric` then argument `formula` should be a full formula with both A and Y. I will just return glm output.")
  } 
  W <- as.matrix(W)
  n <- nrow(W)
  p <- ncol(W)
  if(is.null(sl3_Learner)) {
    if(learning_method == "autoHAL" ) {
      sl3_Learner <- autoML(n,p, parallel, fast_analysis)
    } else if(learning_method == "glmnet" ) {
      sl3_Learner <- Lrnr_glmnet$new()
    } else if(learning_method == "glm" ) {
      sl3_Learner <- Lrnr_glm$new()
    } else if(learning_method == "gam" ) {
      sl3_Learner <- Lrnr_gam$new()
    } else if(learning_method == "mars" ) {
      sl3_Learner <- Lrnr_earth$new()
    } else if (learning_method == "ranger") {
      sl3_Learner <- Lrnr_cv$new(Lrnr_ranger$new())
    } else if(learning_method == "xgboost" ) {
      sl3_Learner <- Stack$new( Lrnr_glmnet$new(), Lrnr_xgboost$new(max_depth =3), Lrnr_xgboost$new(max_depth =4), Lrnr_xgboost$new(max_depth =5), Lrnr_xgboost$new(max_depth =6), Lrnr_xgboost$new(max_depth =7))
      sl3_Learner <- make_learner(Pipeline, Lrnr_cv$new(sl3_Learner), Lrnr_cv_selector$new(loss_squared_error))
    }
    if(cross_fit & learning_method %in% c("glm", "gam", "glmnet", "mars") ) {
      sl3_Learner <- Lrnr_cv$new(sl3_Learner)
    }
  }
   
   
   
  sl3_Lrnr_A <- sl3_Learner
  sl3_Lrnr_Y <- sl3_Learner
  sl3_Lrnr_Y0W <- sl3_Learner
  
  if(learning_method == "glm") {
    sl3_Lrnr_A <- Lrnr_glm$new(formula = glm_formula_A)
    sl3_Lrnr_Y <- Lrnr_glm$new(formula = glm_formula_Y)
    sl3_Lrnr_Y0W <- Lrnr_glm$new(formula = glm_formula_Y0W)
  }
  
  if(estimand == "RR") {
    return(spRR(formula_logRR =  formula, W, A, Y,  sl3_Lrnr_A = sl3_Lrnr_A, sl3_Lrnr_Y = sl3_Lrnr_Y,   weights = weights,  smoothness_order = smoothness_order, max_degree = max_degree, num_knots = num_knots,  fit_control = list(parallel = parallel)))
  }
  if(inference_type == "semiparametric") {
    if(estimand == "CATE") {
      return(spCATE(formula_CATE =  formula, W, A, Y,  sl3_Lrnr_A = sl3_Lrnr_A, sl3_Lrnr_Y = sl3_Lrnr_Y,   weights = weights,  smoothness_order_Y0W = smoothness_order, max_degree_Y0W = max_degree, num_knots_Y0W = num_knots,  fit_control = list(parallel = parallel)))
    }
    if(estimand == "OR") {
      return(spOR(formula_logOR = formula, W, A, Y, weights = weights, W_new = W,  sl3_learner_A = sl3_Lrnr_A, glm_formula_Y0W = glm_formula_Y0W, smoothness_order_Y0W = smoothness_order, max_degree_Y0W = max_degree, num_knots_Y0W = num_knots, reduce_basis = 1e-3, fit_control = list(parallel = parallel), sl3_learner_default = sl3_Lrnr_Y0W    ) )
    }
  } 
  
    
  
}





#'  
autoML <- function(n, p, parallel = F, fast_analysis = F) {
  fit_control <- list(parallel = parallel)
  
  if(p >= 200) {
    lrnr <- Lrnr_glmnet$new()
    return(lrnr)
  }
  
  if(p >= 50) {
    lrnr <-  Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots =  c(5,0), fit_control = fit_control)
    return(lrnr)
  }
  
  if(p >= 25) {
    max_degree <- 1
     
  } else {
    max_degree <- 2
     
  }
   
  if(n<=50) {
    lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots =  c(1,0), fit_control = fit_control)
  }
  else if(n<=100) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(1,1), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 1, num_knots =  c(1,1), fit_control = fit_control)
    }
  } else if(n<=250) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(7,3), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(3,1), fit_control = fit_control)
    }
  } else if(n<=500) {
    if(p<=5) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = 2, num_knots =  c(12,5), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(10,5), fit_control = fit_control)
    }
  } else if(n<=1000) {
    lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(20,10), fit_control = fit_control)
  } else {
    if(fast_analysis) {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(10,5), fit_control = fit_control)
    } else {
      lrnr <- Lrnr_hal9001$new(smoothness_orders = 1, max_degree = max_degree, num_knots =  c(25,15), fit_control = fit_control)
    }
  }
  return(lrnr)
}



