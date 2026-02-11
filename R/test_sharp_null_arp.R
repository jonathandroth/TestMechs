#' @title Hypothesis test for the sharp null
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
#' @param B Bootstrap size, default is zero
#' @param cluster Cluster for bootstrap
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @param alpha Significance level. Default is 0.05
test_sharp_null_arp <- function(df,
                            d,
                            m,
                            y,
                            ordering = NULL,
                            B = 500,
                            cluster = NULL,
                            weight.matrix = "diag",
                            num_Ybins = NULL,
                            alpha = 0.05){

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
  # eta: nuisance par to create target par (= theta_kk TV_kk)

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "zeta" = K,
                   "kappa" = d_y * K, "eta" = K)
  len_x <- sum(par_lengths)

  parameter_types <- c(rep("theta", par_lengths[["theta"]]),
                       rep("delta", par_lengths[["delta"]]),
                       rep("zeta", par_lengths[["zeta"]]),
                       rep("kappa", par_lengths[["kappa"]]),
                       rep("eta", par_lengths[["eta"]]))

  theta_indices <- which(parameter_types == "theta")
  delta_indices <- which(parameter_types == "delta")
  zeta_indices <- which(parameter_types == "zeta")
  kappa_indices <- which(parameter_types == "kappa")
  eta_indices <- which(parameter_types == "eta")



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
    A.shp[k, sum(par_lengths[1:2]) + k] <- 1
    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  }

  #Remove both the extraneous thetas *and* kappa,zeta,eta
  A.shp <- A.shp[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

  #Add shape constraint that all parameters are >= 0 (this is not enforced by ARP)
  A.shp <- rbind(A.shp, diag(NCOL(A.shp)))
  beta.shp <- rep(0, NROW(A.shp))


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

    # Inequalities to equalities
    ## A.obs[2 * K + k, sum(par_lengths[1:2]) + k] <- 1
    A.obs[cbind(2 * K + ((k-1) * d_y + 1):(k * d_y),
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- 1
  }

  #Remove both the extraneous thetas *and* kappa,zeta,eta
  A.obs <- A.obs[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

  #The first 2K rows of A.obs are equality constraints, so we duplicate them with opposite signs to get equalities as inequalities
  A.obs <- rbind(A.obs[1:(2*K),],
                 -A.obs[1:(2*K),],
                 A.obs[(2*K+1):NROW(A.obs), ])

  yvalues <- unique(yvec)
  mvalues <- unique(mvec)
  my_values <- tidyr::expand_grid(m=mvalues,y=yvalues) %>%
                dplyr::arrange(m,y) %>%
                dplyr::select(y,m)

  # Constructing beta_obs
  get_beta.obs <- function(yvec, dvec, mvec) {
    #Get frequencies for all possible values of (y,m) | D=0
    p_ym_0_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                                 .f = ~mean(yvec[dvec == 0] == my_values$y[.x]
                                            & mvec[dvec == 0] == my_values$m[.x]) )

    #Rearrange into a matrix where rows correspond to values of y and cols vals of m
    p_ym_d0 <-
    cbind(my_values, p_ym_0_vec) %>%
      reshape2::dcast(y ~ m, value.var = "p_ym_0_vec") %>%
      dplyr::select(-y) %>%
      as.matrix()


    #Get frequencies for all possible values of (y,m) | D=1
    p_ym_1_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                                 .f = ~mean(yvec[dvec == 1] == my_values$y[.x]
                                            & mvec[dvec == 1] == my_values$m[.x]) )

    #Rearrange into a matrix where rows correspond to values of y and cols vals of m
    p_ym_d1 <-
      cbind(my_values, p_ym_1_vec) %>%
      reshape2::dcast(y ~ m, value.var = "p_ym_1_vec") %>%
      dplyr::select(-y) %>%
      as.matrix()

    # # Matrices (with dimension d_y x K) that store P(y,m | d) for d = 0,1
    # # d = 0
    # p_ym_d0 <- prop.table(table(factor(yvec)[dvec == 0], factor(mvec)[dvec == 0]))
    # # d = 1
    # p_ym_d1 <- prop.table(table(factor(yvec)[dvec == 1], factor(mvec)[dvec == 1]))

    # Matrices that store P(m | d) for d = 0,1
    # d = 0
    p_m_d0 <- colSums(p_ym_d0)
    # d = 1
    p_m_d1 <- colSums(p_ym_d1)

    #Duplicated the first two sets of rows with opposite signs
      # to cast equalities as inequalities
    beta.obs <- c(p_m_d0, p_m_d1,
                  -p_m_d0, -p_m_d1,
                  p_ym_d1 - p_ym_d0)

    return(beta.obs)
  }


  # Bootstrap sample indices
  n <- NROW(df)

  beta.obs_list <- compute_bootstrap_draws_clustered(f =
                                                       function(df,d,y,m,...){get_beta.obs(
                                                                      df[[y]],
                                                                      df[[d]],
                                                                      df[[m]])},
                                                     df = df,
                                                     d = d,
                                                     m = m,
                                                     y = y,
                                                     cluster = cluster,
                                                     numdraws = B,
                                                     return_df = F)



  # Get variance matrix of the beta.obs boostraps
  sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                   beta.obs_list))
  #Get the beta.obs using actual data
  beta.obs <- get_beta.obs(yvec, dvec, mvec)

  #Add the shape constraints as moments
  beta <- c(beta.obs, beta.shp)
  A <- rbind(A.obs, A.shp)
  #Update the covariance matrix to have zero in blocks for the shape constraints
  sigma <- rbind(cbind(sigma.obs, matrix(0,
                                         nrow = NROW(sigma.obs),
                                         ncol = NROW(A.shp) ) ),
                 matrix(0, nrow = NROW(A.shp), ncol = NCOL(sigma.obs) + NROW(A.shp) ))

  #Test if sigma is numerically not psd. If so, add small amount of noise
  min_eig <- base::min(base::eigen(sigma, only.values = TRUE)$values)
  if(min_eig < 0){sigma <- sigma + diag(10*abs(min_eig),
                                        ncol = NROW(sigma),
                                        nrow = NROW(sigma))}


  # Run ARP test
  arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                             y_T = beta,
                                             X_T = A,
                                             sigma = sigma,
                                             alpha = alpha,
                                             hybrid_flag = "ARP"
                                             )

  return(arp)
}
