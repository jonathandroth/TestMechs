#' @title Hypothesis test for the sharp null via Cho and Russell (2023; CR)
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. The inference is via applying CR
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param ordering A list with length equal to the cardinality of the support of
#'   the mediator variable. The name of each element corresponds to a point in
#'   the support, and each element is a vector that collects all m values that
#'   are less than or equal to this point. If ordering = NULL, the standard
#'   ordering is used.
#' @param B Bootstrap size, default is 500
#' @param eps_bar Perturbation parameter used to perturb the
#'   objective/constraints
#' @param alpha Significance level. Default value is .05
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
test_sharp_null_cr <- function(df, d, m, y, ordering = NULL, B = 500,
                               eps_bar = 1e-03, alpha = .05, num_Ybins = NULL){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)


  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
df[[y]] <- yvec
  }
  
  dvec <- df[[d]]
  mvec <- df[[m]]

  # Cardinality of the supports of Y and M
  d_y <- length(unique(yvec))
  K <- length(unique(mvec))

  # Sample size
  n <- nrow(df)

  # Make sure the ordering is defined for all support points
  if (!is.null(ordering) & !all(as.character(unique(mvec)) %in% names(ordering))) {
    stop("Variable ordering does not include all possible values of m!")
  }

  # Define "less than or equal to" function for the given partial ordering
  po_leq <- function(l, k) {
    ifelse(is.null(ordering), l <= k, l %in% ordering[[as.character(k)]])
  }

  # Getting ready to use lpinfer
  # x = (theta, delta, eta)
  # zeta, kappa: nuisance par to convert inequalities to equalities
  # eta: nuisance par to create target par (eta_k = theta_kk TV_kk)

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "eta" = K)
  len_x <- sum(par_lengths)

  # Set theta_lk = 0 for l > k using
  # Matrix that encodes whether l > k
  l_gt_k_mat <- matrix(NA, K, K)
  for (l in 1:K) {
    for (k in 1:K) {
      l_gt_k_mat[l, k] <- !po_leq(l, k)
    }
  }

  # Get hold of the indices
  l_gt_k_inds <- which(l_gt_k_mat)

  # Set shape constraints using A.shp x = beta.shp
  A.shp <- matrix(0, nrow = K, ncol = len_x)

  # Set bound on sum_{l!=k} theta_lk - sum_q delta_qk
  # Also includes eta_kk (= theta_kk TV_kk) as a "nuisance" parameter
  for (k in 1:K) {
    A.shp[k, ((k-1) * K + 1):(k * K)] <- 1
    A.shp[k, (k-1) * K + k] <- 0
    A.shp[k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:2]) + k] <- 1
  }
  A.shp <- A.shp[, -l_gt_k_inds]
  beta.shp <- rep(0, K)


  # Set remaining constraints using A.obs x = beta.obs
  # Define A.obs
  A.obs <- matrix(0,
                  nrow = K + K + (d_y * K),
                  ncol = len_x)

  for (k in 1:K) {
    # Match P(M = k | D = 0)
    A.obs[k, k + K * seq(0, K-1)] <- 1

    # Match P(M = k | D = 1)
    A.obs[K + k, ((k-1) * K + 1):(k * K)] <- 1

    ## # sum theta_{l<k} bound
    ## A.obs[2 * K + k, (k-1) * K + k] <- -1
    ## A.obs[2 * K + k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

    # sup Delta(A) bound
    A.obs[2 * K + ((k-1) * d_y + 1):(k * d_y),
          par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- diag(d_y)

  }

  A.obs <- A.obs[,-l_gt_k_inds]

  # Define beta.obs
  beta.obs <- get_beta.obs(factor(yvec), dvec, factor(mvec))

  # Equalities to two inequalities
  A.obs <- rbind(A.obs, -A.obs[1:(2 * K),])
  beta.obs <- c(beta.obs, -beta.obs[1:(2 * K)])


  # Define target parameter
  A.tgt <- numeric(len_x)
  A.tgt[sum(par_lengths[1:2]) + (1:K)] <- 1
  A.tgt <- A.tgt[-l_gt_k_inds]

  # Update number of parameters
  len_x <- length(A.tgt)

  # Defining model for gurobi::gurobi
  model <- list()

  A <- rbind(A.shp, A.obs)
  rhs <- c(beta.shp, beta.obs)
  lb <- rep(0, len_x)
  ub <- rep(1, len_x)

  # Optimization parameters: suppress output
  params <- list(OutputFlag=0)

  # Combine lower and upper bound into matrix
  model$A <- A
  model$obj <- A.tgt
  model$rhs <- rhs
  model$lb  <- lb
  model$ub  <- ub
  model$sense <- rep('>', length(rhs))

  # Optimize for "l"
  model$modelsense <- 'min'
  min.result <- gurobi::gurobi(model, params)

  # Optimize for "u"
  model$modelsense <- 'max'
  max.result <- gurobi::gurobi(model, params)

  # Draw perturbation
  xi_obj <- runif(len_x) * eps_bar
  xi_rhs <- runif(length(rhs)) * eps_bar
  xi_lb <- runif(len_x) * eps_bar
  xi_ub <- runif(len_x) * eps_bar

  ############################################################################
  # Compute LB-,LB+,UB-,UB+
  ############################################################################

  model$rhs <- rhs - xi_rhs
  model$lb <- lb - xi_lb
  model$ub <- ub + xi_ub

  ############################################################################
  # Compute LB-,UB-
  model$obj<- A.tgt - xi_obj

  # Optimize LB-
  model$modelsense <- 'min'
  min.result.m <- gurobi::gurobi(model, params)

  # Record lower bound
  lbminus <- min.result.m$objval
  if (length(lbminus) == 0) stop("length(lbminus) is zero!")
  # Optimize UB-
  model$modelsense <- 'max'
  max.result.m <- gurobi::gurobi(model, params)

  # Record lower bound
  ubminus <- max.result.m$objval

  #######################################################
  # Compute LB+,UB+
  model$obj <- A.tgt + xi_obj

  # Optimize
  model$modelsense <- 'min'
  min.result.p <- gurobi::gurobi(model, params)

  # Record lower bound
  lbplus <- min.result.p$objval

  # Optimize UB-
  model$modelsense <- 'max'
  max.result.p <- gurobi::gurobi(model, params)

  # Rrecord upper bound
  ubplus<- max.result.p$objval

  ############################################################################
  # Begin bootstrap procedure
  boot_lbminus <- rep(NA, B)
  boot_lbplus <- rep(NA, B)
  boot_ubminus <- rep(NA, B)
  boot_ubplus <- rep(NA, B)

  for (b in 1:B) {

    # Bootstrap indices
    b_ind <- sample.int(n, n, replace = T)

    # Update beta.obs
    beta.obs_b <- get_beta.obs(factor(yvec)[b_ind],
                               dvec[b_ind],
                               factor(mvec)[b_ind])


    beta.obs_b <- c(beta.obs_b, -beta.obs_b[1:(2 * K)])
    rhs <- c(beta.shp, beta.obs_b)

    # Update model
    model$rhs <- rhs - xi_rhs

    ############################################################################
    # Compute LB-,UB-
    model$obj<- A.tgt - xi_obj

    # Optimize LB-
    model$modelsense <- 'min'
    bmin.result.m <- gurobi::gurobi(model, params)

    # Record lower bound
    boot_lbminus[b] <- bmin.result.m$objval

    # Optimize UB-
    model$modelsense <- 'max'
    bmax.result.m <- gurobi::gurobi(model, params)

    # Record upper bound
    boot_ubminus[b] <- bmax.result.m$objval

    ###################################################
    # Objective function for LB+ and UB+
    model$obj<- A.tgt + xi_obj

    # Optimize LB+
    model$modelsense <- 'min'
    bmin.result.p <- gurobi::gurobi(model, params)

    # Record lower bound
    boot_lbplus[b] <- bmin.result.p$objval

    # Optimize UB+
    model$modelsense <- 'max'
    bmax.result.p <- gurobi::gurobi(model, params)

    # Record upper bound
    boot_ubplus[b] <- bmax.result.p$objval
  }

  # Compute indicator Dn
  bn <- 1/sqrt(log(n))
  Delta <- max(ubplus, ubminus) - min(lbminus, lbplus)
  Dn <- (Delta > bn) + 0

  # Calculate kappa
  kappa <- (1 - alpha) * Dn + (1 - alpha/2) * (1 - Dn)

  # Bootstrap quantities
  lbminus_q <- sqrt(n) * (boot_lbminus - lbminus)
  lbplus_q <- sqrt(n) * (boot_lbplus - lbplus)
  ubminus_q <- -sqrt(n) * (boot_ubminus - ubminus)
  ubplus_q <- -sqrt(n) * (boot_ubplus - ubplus)

  # Select quantile according to kappa
  psi_k_lb_minus <- quantile(lbminus_q, kappa)
  psi_k_lb_plus <- quantile(lbplus_q, kappa)
  psi_k_ub_minus <- quantile(ubminus_q, kappa)
  psi_k_ub_plus <- quantile(ubplus_q, kappa)

  #compute confidence set for alpha=0.05
  CSlb <- min(lbminus, lbplus) - (1/sqrt(n)) * max(psi_k_lb_minus, psi_k_lb_plus)
  CSub <- max(ubminus, ubplus) + (1/sqrt(n)) * max(psi_k_ub_minus, psi_k_ub_plus)

  return(list(CI = c(CSlb, CSub), reject = (0 < CSlb)))
}
