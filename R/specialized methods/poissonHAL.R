
##### Need to install hal9001: devtools::install_github("tlverse/hal9001", ref = "devel")
library(hal9001)
library(data.table)
library(glmnet)
############################
#### EXAMPLE ################
############################
# n <- 1000
# W <-  runif(n,-1, 1)
# W1 <-   runif(n,-1, 1)
# A <-    rbinom(n,size=1, 1/2)
# Ttilde <- rweibull(n, shape = 2, scale = 1/exp(A - +W1+ W))
# Delta <- rep(1, n)
# X <- cbind(W,A,W1)
#
#
# haz <- dweibull(Ttilde, shape = 2, scale = 1/exp(A - W)) / (1-pweibull(Ttilde, shape = 2, scale = 1/exp(A - W)))
#
#
#
# fit <- fit_hal_hazard(X, Ttilde, Delta, num_bins_X = 10, num_bins_t = 20, formula_hal = "~ . + h(t,.) + h(A,.) + h(A,t,.)",   penalize_baseline_time = F, parallel = TRUE, ncores = 4, max_num_groups =  500  )
#
#
# print(as.data.table(fit$X_poisson))
# print(length(fit$basis_list))
# print(length(fit$basis_list_t))
# print(length(fit$basis_list_t_X))
#
# drop <- Ttilde > quantile(Ttilde, 0.9)
# preds <- as.vector(predict(fit, X, Ttilde, each_t = F))
#
# plot(Ttilde[-drop], preds[-drop])
# plot(Ttilde[-drop], haz[-drop])


#' @param X A matrix of baseline variables
#' @param Ttilde A vector of event times. Ttilde = min(T,C)
#' @param Delta A binary vector where Delta=1 is the event type of interest. Delta = T <= C
#' @param num_bins_t Total numbers of bins to use for the time component.
#' @param num_bins_X A vector of interaction degree specific Maximum number of bins to use for each covariate when making the time-independent HAL basis.
#' @param max_degree_XT Maximum interaction degree of basis functions. max_degree_XT = 1 means no interactions between time and X.
#' max_degree_XT = 2 generates up to a 2-tensor product of X basis functions and time basis functions.
#' max_degree_XT = 3 generates up to a 3-tensor product containing 2 main-term X basis functions and 1 main-term time basis functions.
#' @param time_grid Grid of times that specifies time bins to use in fitting. By default, a quantile-equal-spaced time-grid is made with num_bins_t knots.
#' @param penalize_baseline_time Whether to penalize the baseline time basis functions.
#' @param formula_hal See the formula argument from the package Github/tlverse/HAL9001. Note this must be a string.
#' @export
#'
#'Useful formulas: formula_hal =
#' "~ ."
#' "~ .^2"
#' "~ . + h(t,.) "
#' "~ . + h(t,.) + h(A,.) + h(A,t,.)" (E.g. A is binary treatment)
#' "~ .^2 + h(t, ., .)"
fit_hal_hazard = function(X, Ttilde, Delta, formula_hal = "~ . + h(t,.) ", num_bins_t = 20, num_bins_X = 10,   max_degree_XT = 2,   time_grid = NULL, weight_mat_by_time_bin = NULL,   penalize_baseline_time = FALSE,  dont_fit = FALSE,   ...) {
  formula_hal <- as.character(formula_hal)
  order = 0
  include_time_interaction = TRUE


  X <- as.matrix(X)
  ################################################################################################
  ######## Generate time grid for binning if not provided ########################################
  ################################################################################################
  if(is.null(time_grid)) {
    # Grid based on event times only????????????
    time_grid = sort(unique(quantile(Ttilde[Delta==1] , seq(0,1, length = num_bins_t))))
    #time_grid <- union( min(Ttilde ), time_grid, max(Ttilde ))
    #time_grid <- time_grid[-length(time_grid)]
  }
  time_grid[which.min(time_grid)] <- min(Ttilde )
  time_grid[which.max(time_grid)] <- max(Ttilde )


  #A <- as.vector(A)
  Ttilde <- as.vector(Ttilde)
  Delta <- as.vector(Delta)
  covariates <- colnames(X)
  if(is.null(covariates)) {
    covariates <- paste0("X",seq_len(ncol(X)))
    colnames(X) <- covariates
  }
  data <- data.frame(X,Ttilde = Ttilde, Delta = Delta)


  X <- as.matrix(data[,covariates, drop = F])

  ################################################################################################
  ######## Generate basis functions based on discretized covariates/time ##########################
  ################################################################################################
  X_discr <- hal9001:::quantizer(X, max(num_bins_X))
  colnames(X_discr) <- covariates
  Ttilde_discr <- as.matrix(time_grid[findInterval(Ttilde, time_grid, all.inside = T )])
  colnames(Ttilde_discr) <- "Ttilde"

  XT_discr <- cbind(X_discr, Ttilde_discr)
  print(colnames(XT_discr))

  ###### Generate minimal basis list that encodes the covariate-space partitioning (We only need main-term basis functions)
  ###### This is just for computational benefits. basis_list_key could be replaced with basis_list_X
  basis_list_key <-  hal9001::enumerate_basis(X_discr, max_degree =1, smoothness_orders = 0,  num_knots = nrow(X_discr)+100)
  basis_list <- basis_list_key
  x_basis <- hal9001::make_design_matrix(X, basis_list_key)
  tx_basis <- Matrix::t(x_basis)
  ######### Generate list/map that maps a unique group index to group members. By group I mean an element/bin of the covariate-space partition.
  row_copy_map <- hal9001::make_copy_map(tx_basis)
  ngroups <- length(row_copy_map)
  n_per_group <- sapply(row_copy_map, length)
  ##### Find unique HAL basis realizations and subset to those only.
  x_basis_key <- x_basis[as.numeric(names(row_copy_map)), ,drop = F]

  x_basis_key <- x_basis[as.numeric(names(row_copy_map)), ,drop = F]
  X_key <- X[as.numeric(names(row_copy_map)),, drop = F]


  if(is.null(weight_mat_by_time_bin)) {
    weight_mat_by_time_bin <- matrix(1, nrow = nrow(X), ncol = length(time_grid)-1 )
  }
  ###### Generate basis functions from binned covariates values
  ###### Generate X-specific basis functions (time independent)
  ##### Generate basis functions that also include interactions between X and T (duplicates with basis_list_X will be removed later)
  if(!is.null(formula_hal)){
    basis_list_X_T <- NULL
    tryCatch({
      basis_list_X_T <<- formula_hal(formula_hal, X = XT_discr, smoothness_orders = order, num_knots = 10000)$basis_list
    },
    error = function(e) {
      XT_discr <- XT_discr
      colnames(XT_discr) <- c(covariates, "t")
      basis_list_X_T <<- formula_hal(formula_hal, X = XT_discr, smoothness_orders = order, num_knots = 10000)$basis_list
    })
    #basis_list_X_T <- formula_hal(formula_hal, X = XT_discr, smoothness_orders = 0 num_knots = 10000)$basis_list
    basis_list_X <- basis_list_X_T
  } else {
    #basis_list_X <- hal9001::enumerate_basis(X_discr, max_degree =max_degree_XT , smoothness_orders = 0, num_knots = 1000)
    basis_list_X_T <- hal9001::enumerate_basis(XT_discr, max_degree =max_degree_XT, smoothness_orders = order, num_knots = 1000)
    basis_list_X <- basis_list_X_T
  }


  time_grid_start <- time_grid[-length(time_grid)]
  ###### Convert weight matrix to reduced dimension by averaging
  weights  <- unlist(lapply(seq_along(time_grid_start) , function(t_index) {



    out <- lapply(row_copy_map, function(grp) {
      mean(weight_mat_by_time_bin[,t_index])
    })
    names(out) <- paste0("grp_",seq_along(row_copy_map))
    return(out)
  }))

  ###### Compute total risk time within each group + time bin

  time_at_risk <- lapply(seq_along(time_grid_start) , function(t_index) {
    tr <- time_grid[t_index]
    tr1 <- time_grid[t_index +1]


    out <- lapply(row_copy_map, function(grp) {
      sum(( Ttilde[grp] >= tr)*(pmin( Ttilde[grp], tr1) - tr))
    })
    names(out) <- paste0("grp_",seq_along(row_copy_map))
    return(out)
  })
  names(time_at_risk) <- paste0("t_", time_grid_start)
  time_at_risk <- unlist(time_at_risk)


  ###### Compute total number of events within each group + time bin
  num_events <- lapply(seq_along(time_grid_start) , function(t_index) {
    tr <- time_grid[t_index]
    tr1 <- time_grid[t_index +1]
    if(t_index +1 == length(time_grid)) {
      tr1 <- Inf
    }

    out <- lapply(row_copy_map, function(grp) {
      sum(Delta[grp]*(Ttilde[grp] >= tr)*(Ttilde[grp] < tr1))
    })
    names(out) <- paste0("grp_", seq_along(row_copy_map))
    return(out)
  })
  names(num_events) <- paste0("t_", time_grid_start)
  num_events <- unlist(num_events)

  if(any(is.na(num_events)) || any(is.na(time_at_risk))) {
    warning("Different NA rows for num_events and time_at_risk")
  }



  ###### Generate reduced datasets with rows corresponding to partition groups + time
  X_poisson <- as.data.table(X_key)
  colnames(X_poisson) <- covariates
  print(time_grid_start)
  X_poisson <- rbindlist(lapply( (time_grid_start), function(t) {
    X_poisson <- copy(X_poisson)
    X_poisson$t <- t
    return(X_poisson)
  }))

  set(X_poisson,,"num_events", num_events)
  set(X_poisson,,"time_at_risk", time_at_risk)
  set(X_poisson ,, "weights", weights)
  X_poisson <- as.matrix(X_poisson)



  # Generate baseline time-only basis function list for fitting HAL

  t_index <- which(colnames(X_poisson)=="t")
  basis_list_t <- hal9001::enumerate_basis(as.matrix(time_grid), max_degree =1, num_knots = length(time_grid)+10, smoothness_orders = order )
  basis_list_t <- lapply(basis_list_t, function(basis) {
    basis$cols <- t_index
    basis
  })



  # Remove empty bins. This code might be dangerous.
  remove_rows <- which(X_poisson[,"time_at_risk"]==0)
  print(remove_rows)
  print(as.data.table(X_poisson))
  if(length(remove_rows) > 0) {
    X_poisson <- X_poisson[-remove_rows,]
  }
  #### Generate unpenalized time-dependent baseline basis matrix.
  x_basis_t_unpenal <- hal9001::make_design_matrix(as.matrix(X_poisson), basis_list_t)

  # Set up poisson HAL
  offset <- log(X_poisson[,"time_at_risk"])
  Y <- X_poisson[,"num_events"]
  X <- X_poisson[,c(covariates,"t")]
  X_unpenalized <-  (x_basis_t_unpenal)
  # Remove constant baseline basis functions or else unpenalized might error
  drop_times <- which(apply(X_unpenalized,2, function(x){
    length(unique(x))
  })==1)
  print("dropa")
  print(drop_times)


  ##### Merge the basis lists for X, X_T, and T
  full_basis_list <- union( basis_list_X, basis_list_X_T)
  full_basis_list <- setdiff(full_basis_list, basis_list_t[-drop_times])
  if(penalize_baseline_time) {
    penalty <- c(rep(1, length(full_basis_list)),  rep(0.5, length(basis_list_t[-drop_times])))
  } else {
    penalty <- c(rep(1, length(full_basis_list)),   rep(0.01,  length(basis_list_t[-drop_times])))
  }
  full_basis_list <- c(full_basis_list, basis_list_t[-drop_times])

  # MAKE HAL design matrix for training
  x_basis <- hal9001::make_design_matrix(X, full_basis_list)
  print("x_basis dimension:")
  print(dim(x_basis))
  print("fitting glmnet")


  if(dont_fit) {
    hal_fit <- list(coefs = NULL, basis_list = full_basis_list, X_poisson = X_poisson, covariates = covariates, X_key = as.data.table(X_key), X_basis_key = x_basis_key, X_basis_list_key = basis_list_key, basis_list_t = basis_list_t, basis_list_X =basis_list_X, basis_list_t_X = basis_list_X_T, penalty = penalty, time_grid = time_grid )
    class(hal_fit) <- "hal_fit_hazard"
    return(hal_fit)
  }
  #### FIT HAL
  # if(parallel) {
  #   doMC::registerDoMC(ncores)
  # }
  # if(no_penalization) {
  #   fit <- speedglm::speedglm.wfit(Y, x_basis, offset = offset, family = poisson())
  #   coefs <- as.vector(coef(fit ))
  # } else {
  fit <- glmnet(x_basis, Y, offset = offset, family = "poisson", penalty.factor  = penalty, parallel  = parallel, weights =weights ,...)
  coefs <- as.vector(coef(fit,"lambda.min"))
# }

  fit_info <- list(glmnet.fit = fit, x_basis_train = x_basis, X_poisson = X_poisson,  X_key = as.data.table(X_key), X_basis_key = x_basis_key, X_basis_list_key = basis_list_key, basis_list_t = basis_list_t, basis_list_X =basis_list_X, basis_list_t_X = basis_list_X_T, penalty = penalty)
  hal_fit <- list(  coefs = coefs, basis_list = full_basis_list,  covariates = covariates,  time_grid = time_grid, fit_info = fit_info )
  class(hal_fit) <- "hal_fit_hazard"
#
#
#   t_grid <- unlist(sapply(1:length(time_grid), function(i) {
#     l <- signif( time_grid[i],3)
#     if(i == length( time_grid)) {
#       return(paste0("(", l,", ", Inf, ")"))
#     }
#     r <- signif( time_grid[i+1],3)
#     if(i==1) {
#       return(c(paste0("(", -Inf,", ", l, "]"), paste0("[", l,", ", r, ")")))
#     }
#     if(i == length( time_grid)-1) {
#       return(paste0("[", l,", ", r, "]"))
#     } else{
#       return(paste0("[", l,", ", r, ")"))
#     }
#
#   }))
#
#
#
#   lambda_t <- (predict(hal_fit, X,  c(-Inf, time_grid), each_t = T))
#   t_diffs <- c(0,diff(time_grid), 0)
#   colnames(lambda_t) <- t_grid
#
#   St <- as.data.table(t(apply(lambda_t, 1, function(haz) {
#     exp(-cumsum(haz * t_diffs))
#   })))
#   St[, ncol(St)] <- 0
#   colnames(St) <- paste0(as.character(c(0, round(time_grid,3))), ">=")
#
#   hal_fit$haz <- lambda_t
#   hal_fit$surv <- St


  return(hal_fit)

}






get_group_index <- function(x, x_basis_key_transpose) {
  index <- which.min(colSums(abs(x - x_basis_key_transpose)))[1]

  return(index)
}


#' @export
predict.hal_fit_hazard <- function(fit, X_new, t_new, each_t = FALSE) {
  if(!is.matrix(X_new)){
    X_new <- as.matrix(X_new)
  }

  t_new <- as.vector(t_new)
  try({X_new <- X_new[, fit$covariates, drop = F]}, silent = TRUE)

  #X_new_basis <- hal9001::make_design_matrix(X_new ,fit$X_basis_list_key)
  #new_grps <- as.vector(apply(X_new_basis, 1, get_group_index, x_basis_key = as.matrix(t(fit$X_basis_key))))
  #fit$X_key <- as.data.table(fit$X_key)
  if(!each_t) {
    #new_X_poisson <- fit$X_key[new_grps, ]
    new_X_poisson <- as.data.table(X_new)
    new_X_poisson$t <- t_new


  } else {
    # stacked_grps <- rep(new_grps, length(t_new))
    # stacked_t <- rep(t_new, each = length(new_grps))
    #new_X_poisson <- fit$X_key[stacked_grps, ]
    stacked_t <- rep(t_new, each = nrow(X_new))
    new_X_poisson <- as.data.table(X_new)
    new_X_poisson <- new_X_poisson[rep(1:nrow(new_X_poisson) ,  length(t_new))]
    new_X_poisson$t <- stacked_t

  }



  #new_X_poisson$t <- t_new
  coefs <- fit$coefs
  x_basis_new <- hal9001::make_design_matrix(as.matrix(new_X_poisson), fit$basis_list, p_reserve=1)
  link <- as.vector(Matrix::tcrossprod(
    x = x_basis_new,
    y = coefs[-1]
  ) + coefs[1])
  predictions <- exp(link)
  if(each_t) {
    predictions <- matrix(predictions, nrow = nrow(X_new), byrow = F)
    colnames(predictions) <- t_new
    t_mat <- matrix(rep(t_new, each = nrow(predictions)), nrow = nrow(predictions))
    keep <- (t_mat < Inf) * (t_mat >= min(fit$time_grid))
    predictions <-predictions*keep

  }
  return(predictions)
}

