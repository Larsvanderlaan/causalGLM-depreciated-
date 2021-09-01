#'
#' @importFrom R6 R6Class
#' @importFrom origami folds2foldvec
#' @importFrom stats predict quasibinomial
#' @importFrom sl3 Lrnr_base
#' @import data.table
#' @import hal9001
#' @import glmnet
#' @export
Lrnr_hal9001_custom <- R6Class(
  classname = "Lrnr_hal9001_custom", inherit = Lrnr_base,
  portable = TRUE, class = TRUE,
  public = list(
    initialize = function(max_degree = 3,
                          fit_type = "glmnet",
                          n_folds = 10,
                          use_min = TRUE,
                          reduce_basis = NULL,
                          return_lasso = TRUE,
                          return_x_basis = FALSE,
                          basis_list = NULL,
                          cv_select = TRUE,
                          formula_hal = NULL,
                          fit_control = list(),
                          ...) {
      params <- sl3:::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    .properties = c("continuous", "binomial", "weights", "ids"),

    .train = function(task) {
      args <- self$params

      outcome_type <- self$get_outcome_type(task)

      if (is.null(args$family)) {
        args$family  <- outcome_type$glm_family()
      }


      if(args$family=="binomial" && !all(task$Y %in% c(0,1))) {
        args$family <- binomial()
        args$Y <- task$Y
      } else {
        args$Y <- outcome_type$format(task$Y)
      }



      args$X <- as.matrix(task$X)

      args$yolo <- FALSE
      args$formula <- args$formula_hal
      #args$fit_control <- list()
      if (task$has_node("weights")) {
        args$fit_control$weights <- task$weights
      }

      if (task$has_node("offset")) {
        args$offset <- task$offset
      }

      if (task$has_node("id")) {
        args$id <- task$id
      }

      # pass in formals of glmnet versus cv.glmnet based on cv_select
      if (args$cv_select) {
        glmnet_other_valid <- union(
          names(formals(glmnet::cv.glmnet)),
          names(formals(glmnet::glmnet))
        )
      } else {
        glmnet_other_valid <- names(formals(glmnet::glmnet))
      }

      # fit HAL, allowing glmnet-fitting arguments
      fit_object <- sl3:::call_with_args(
        hal9001::fit_hal, args,
        other_valid = glmnet_other_valid
      )
      return(fit_object)
    },
    .predict = function(task = NULL) {
      predictions <- predict(self$fit_object, new_data = as.matrix(task$X))
      if (!is.na(safe_dim(predictions)[2])) {
        p <- ncol(predictions)
        colnames(predictions) <- sprintf("lambda_%0.3e", self$params$lambda)
      }
      return(predictions)
    },
    .required_packages = c("hal9001", "glmnet")
  )
)
