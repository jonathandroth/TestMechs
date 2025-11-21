#' @title Hypothesis test for the sharp null via Fang et al. (2023; FSST)
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. The inference is via applying FSST
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
#' @param B Bootstrap size, default is zero
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param alpha Significance level. Default value is .05.
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
test_sharp_null_fsst <- function(df, d, m, y, ordering = NULL, B = 500,
                                 weight.matrix = "diag", alpha = .05,
                                 num_Ybins = NULL, lambda = NA){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)


  yvec <- df[[y]]
  n <- length(yvec)
  
  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }


  dvec <- df[[d]]
  mvec <- df[[m]]

  # Cardinality of the supports of Y and M
  d_y <- length(unique(yvec))
  K <- length(unique(mvec))

  # Make sure the ordering is defined for all support points
  if (!is.null(ordering) & !all(as.character(unique(mvec)) %in% names(ordering))) {
    stop("Variable ordering does not include all possible values of m!")
  }

  # Define "less than or equal to" function for the given partial ordering
  po_leq <- function(l, k) {
    ifelse(is.null(ordering), l <= k, l %in% ordering[[as.character(k)]])
  }

  # Getting ready to use lpinfer
  # x = (theta, delta, zeta, kappa, eta)
  # zeta, kappa: nuisance par to convert inequalities to equalities
  # eta: nuisance par to create target par (eta_k = theta_kk TV_kk)

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "zeta" = K,
                   "kappa" = d_y * K, "eta" = K)
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

    # Inequalities to equalities
    A.shp[k, sum(par_lengths[1:2]) + k] <- -1
    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  }
  A.shp <- A.shp[, -l_gt_k_inds]
  beta.shp <- rep(0, K)


  # Set remaining constraints using A.obs x = beta.obs

  # Bootstrap sample indices
  boot_mat <- matrix(sample(1:nrow(df), B * nrow(df), replace = T), nrow = B)
  # Get all beta_obs to pass to lpinfer
  beta.obs_list <- lapply(1:B,
                          function(b) get_beta.obs(factor(yvec)[boot_mat[b,]],
                                                   dvec[boot_mat[b,]],
                                                   factor(mvec)[boot_mat[b,]])
                          )

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

    # Inequalities to equalities
    ## A.obs[2 * K + k, sum(par_lengths[1:2]) + k] <- 1
    A.obs[cbind(2 * K + ((k-1) * d_y + 1):(k * d_y),
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- -1
  }

  A.obs <- A.obs[,-l_gt_k_inds]

  # Define target parameter
  A.tgt <- numeric(len_x)
  A.tgt[sum(par_lengths[1:4]) + (1:K)] <- 1
  A.tgt <- A.tgt[-l_gt_k_inds]

  beta.obs <- get_beta.obs(factor(yvec), dvec, factor(mvec))
  # Run FSST
  lpm <- lpinfer::lpmodel(A.obs = A.obs,
                          A.shp = A.shp,
                          A.tgt = A.tgt,
                          beta.obs = c(list(beta.obs),
                                       beta.obs_list),
                          beta.shp = beta.shp)

  if (!is.na(lambda)) {
    lambda <- 1/sqrt(log(max(length(beta.obs), exp(1))) * log(max(exp(1), log(max(exp(1), n)))))
  } else {
    lambda <- NA
  }
  
  fsst_result <- lpinfer::fsst(df, lpmodel = lpm, beta.tgt = 0, R = B,
                               weight.matrix = weight.matrix, lambda = lambda)

  return(list(result = fsst_result, reject = (fsst_result$pval[,2] < alpha)))
}

#' @title Hypothesis test for the sharp null with alternative LP formulation
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. The inference is via applying Fang et al. (2023; FSST)
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
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param alpha Significance level. Default value is .05.
test_sharp_null_alt <- function(df,
                                d,
                                m,
                                y,
                                ordering = NULL,
                                B = 500,
                                weight.matrix = "diag",
                                alpha = .05){

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  # Cardinality of the supports of Y and M
  d_y <- length(unique(yvec))
  K <- length(unique(mvec))

  # Make sure the ordering is defined for all support points
  if (!is.null(ordering) & !all(as.character(unique(mvec)) %in% names(ordering))) {
    stop("Variable ordering does not include all possible values of m!")
  }

  # Define "less than or equal to" function for the given partial ordering
  po_leq <- function(l, k) {
    ifelse(is.null(ordering), l <= k, l %in% ordering[[as.character(k)]])
  }

  # Getting ready to use lpinfer
  # x = (theta, delta, zeta, kappa, eta)
  # zeta, kappa: nuisance par to convert inequalities to equalities
  # eta: nuisance par to create target par (eta_k = theta_kk TV_kk)

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "zeta" = K,
                   "kappa" = d_y * K, "eta" = K)
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
  A.shp <- matrix(0, nrow = length(l_gt_k_inds), ncol = len_x)

  A.shp[cbind(1:length(l_gt_k_inds), l_gt_k_inds)] <- 1
  # Set bound on sum_{l!=k} theta_lk - sum_q delta_qk
  # Also includes eta_kk (= theta_kk TV_kk) as a "nuisance" parameter
  ## for (k in 1:K) {
  ##   A.shp[k, ((k-1) * K + 1):(k * K)] <- 1
  ##   A.shp[k, (k-1) * K + k] <- 0
  ##   A.shp[k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

  ##   # Inequalities to equalities
  ##   A.shp[k, sum(par_lengths[1:2]) + k] <- -1
  ##   # TV_kk "nuisance" parameters
  ##   A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  ## }
  ## A.shp <- A.shp[, -l_gt_k_inds]
  beta.shp <- rep(0, length(l_gt_k_inds))


  # Set remaining constraints using A.obs x = beta.obs

  # Bootstrap sample indices
  boot_mat <- matrix(sample(1:nrow(df), B * nrow(df), replace = T), nrow = B)
  # Get all beta_obs to pass to lpinfer
  beta.obs_list <- lapply(1:B,
                          function(b) get_beta.obs_alt(factor(yvec)[boot_mat[b,]],
                                                       dvec[boot_mat[b,]],
                                                       factor(mvec)[boot_mat[b,]])
                          )

  # Define A.obs
  A.obs <- matrix(0,
                  nrow = K + K + K + (d_y * K),
                  ncol = len_x)

  for (k in 1:K) {
    # Match P(M = k | D = 0)
    A.obs[k, k + K * seq(0, K-1)] <- 1

    # Match P(M = k | D = 1)
    A.obs[K + k, ((k-1) * K + 1):(k * K)] <- 1

    ## # sum theta_{l<k} bound
    A.obs[2 * K + k, (k-1) * K + k] <- -1
    A.obs[2 * K + k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

    # sup Delta(A) bound
    A.obs[3 * K + ((k-1) * d_y + 1):(k * d_y),
          par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- diag(d_y)

    # Inequalities to equalities
    A.obs[2 * K + k, sum(par_lengths[1:2]) + k] <- -1
    A.obs[cbind(3 * K + ((k-1) * d_y + 1):(k * d_y),
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- -1
  }

  ## A.obs <- A.obs[,-l_gt_k_inds]

  # Define target parameter
  A.tgt <- numeric(len_x)
  A.tgt[sum(par_lengths[1:4]) + (1:K)] <- 1
  ## A.tgt <- A.tgt[-l_gt_k_inds]

  # Run FSST
  lpm <- lpinfer::lpmodel(A.obs = A.obs,
                          A.shp = A.shp,
                          A.tgt = A.tgt,
                          beta.obs = beta.obs_list,
                          beta.shp = beta.shp)

  fsst_result <- lpinfer::fsst(df, lpmodel = lpm, beta.tgt = 0, R = B-1,
                               weight.matrix = weight.matrix)

  return(list(result = fsst_result, reject = (fsst_result$pval < alpha)))
}


#' @title Function to construct beta.obs to pass to lpinfer
##' @param yvec Vector of outcome variable
##' @param dvec Vector of treatment variable
##' @param mvec Vector of mediator variable
get_beta.obs <- function(yvec, dvec, mvec) {

  # Matrices (with dimension d_y x K) that store P(y,m | d) for d = 0,1
  # d = 0
  p_ym_d0 <- prop.table(table(yvec[dvec == 0], mvec[dvec == 0]))
  # d = 1
  p_ym_d1 <- prop.table(table(yvec[dvec == 1], mvec[dvec == 1]))

  # Matrices that store P(m | d) for d = 0,1
  # d = 0
  p_m_d0 <- colSums(p_ym_d0)
  # d = 1
  p_m_d1 <- colSums(p_ym_d1)

  beta.obs <- c(p_m_d0, p_m_d1, p_ym_d1 - p_ym_d0)

  return(beta.obs)
}

#' @title Function to construct beta.obs to pass to lpinfer under an alternative
#'   LP formulation
##' @param yvec Vector of outcome variable
##' @param dvec Vector of treatment variable
##' @param mvec Vector of mediator variable
get_beta.obs_alt <- function(yvec, dvec, mvec) {

  # Matrices (with dimension d_y x K) that store P(y,m | d) for d = 0,1
  # d = 0
  p_ym_d0 <- prop.table(table(yvec[dvec == 0], mvec[dvec == 0]))
  # d = 1
  p_ym_d1 <- prop.table(table(yvec[dvec == 1], mvec[dvec == 1]))

  # Matrices that store P(m | d) for d = 0,1
  # d = 0
  p_m_d0 <- colSums(p_ym_d0)
  # d = 1
  p_m_d1 <- colSums(p_ym_d1)

  beta.obs <- c(p_m_d0, p_m_d1, -p_m_d1, p_ym_d1 - p_ym_d0)

  return(beta.obs)
}
