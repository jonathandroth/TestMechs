#' @title Hypothesis test for the sharp null
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. The mediator is assuemd to be binary. 
#'   The inference is via Cox and Shi (2023) OR Andrews, Roth, and Pakes OR FSST.
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param method A string variable indicating which method: CS, ARP, FSST
#' @param ordering A list with length equal to the cardinality of the support of
#'   the mediator variable. The name of each element corresponds to a point in
#'   the support, and each element is a vector that collects all m values that
#'   are less than or equal to this point. If ordering = NULL, the standard
#'   ordering is used.
#' @param B Bootstrap size, default is zero
#' @param cluster Cluster for bootstrap
#' @param reg_formula Regression formula for the non-experimental setting
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#'  @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#'  @param fix_n1 Whether the number of treated units (or clusters) should be fixed in the bootstrap
#'  @param lambda (For FSST only) A string variable, either dd or ndd, standing for data-driven or non-data driven respectively.
#'  @param analytic_variance (For CS or ARP only) A flag indicating whether to use analytic variance.
#'  @param refinement (For CS only, optional) If TRUE, use the refined Cox & Shi test (rCC rather than CC). Default is FALSE.
#'  @param print_both_var (For CS only, Optional) If TRUE, print sigma from both actual and bootstrapped data.
#' @export

test_sharp_null_binary_m <- function(df,
                                     d,
                                     m,
                                     y,
                                     method = NULL,
                                     ordering = NULL,
                                     B = 500,
                                     cluster = NULL,
                                     reg_formula = NULL,
                                     weight.matrix = "diag",
                                     ats_only = F,
                                     alpha = 0.05,
                                     kappa = alpha/10,
                                     use_hybrid = T,
                                     num_Ybins = NULL,
                                     fix_n1 = T,
                                     lambda = "dd",    #fsst arg
                                     analytic_variance = FALSE,    # arp cs arg
                                     refinement = FALSE,    #cs arg
                                     print_both_var = FALSE    #cs arg
                                     ){
  
  ## Remove missing
  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)
  
  ## Pre-processing
  yvec <- df[[y]]
  n <- length(yvec)
  
  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  } else {
    continuous_y_flag <- n / length(unique(yvec)) <= 30
    if (continuous_y_flag) {
      message("Y variable might be continuous. Discretize it by specifying num_Ybins. Default num_Ybins = 5 is used now.")
      num_Ybins <- 5
      yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
      df[[y]] <- yvec
    }
  }
  
  if (is.null(cluster)) {
    clustervec <- 1:length(yvec)
  } else {
    clustervec <- df[[cluster]]
  }
  dvec <- df[[d]]
  mvec <- df[[m]]
  
  yvalues <- sort(unique(yvec))
  mvalues <- unique(mvec)
  my_values <- purrr::cross_df(list(m=mvalues,y=yvalues)) %>%
    dplyr::arrange(m,y) %>%
    dplyr::select(y,m)
  
  
  ## Define function for beta.obs
  get_beta.obs <- function(yvec, dvec, mvec, df, d, yvalues, reg_formula) {
    if (!is.null(reg_formula)) {
      # Use shared regression helpers from test_sharp_null.R
      p_ym_0_vec <- compute_regression_probs(df, yvec, mvec, my_values, d, reg_formula, control_transform)
      p_ym_1_vec <- compute_regression_probs(df, yvec, mvec, my_values, d, reg_formula, treated_transform)
    }

    # Get partial density for Y,M = 1 | D = d
    if (is.null(reg_formula)) {
      p_y1_1 <- purrr::map_dbl(.x = 1:length(yvalues),
                               .f = ~mean(yvec[dvec == 1] == yvalues[.x]
                                          & mvec[dvec == 1] == 1 ))
      p_y1_0 <- purrr::map_dbl(.x = 1:length(yvalues),
                               .f = ~mean(yvec[dvec == 0] == yvalues[.x]
                                          & mvec[dvec == 0] == 1 ))
    } else {
      # Rows with M = 1 map to the p_y1_1 / p_y1_0 probabilities
      p_y1_1 <- p_ym_1_vec[my_values$m == 1]
      p_y1_0 <- p_ym_0_vec[my_values$m == 1]
    }

    if (!ats_only) {
      # Get partial density for Y,M = 0 | D = d
      if (is.null(reg_formula)) {
        p_y0_1 <- purrr::map_dbl(.x = 1:length(yvalues),
                                 .f = ~mean(yvec[dvec == 1] == yvalues[.x]
                                            & mvec[dvec == 1] == 0 ))
        p_y0_0 <- purrr::map_dbl(.x = 1:length(yvalues),
                                 .f = ~mean(yvec[dvec == 0] == yvalues[.x]
                                            & mvec[dvec == 0] == 0 ))
      } else {
        # Rows with M = 0 map to the p_y0_1 / p_y0_0 probabilities
        p_y0_1 <- p_ym_1_vec[my_values$m == 0]
        p_y0_0 <- p_ym_0_vec[my_values$m == 0]
      }
    }
    #We return differences in partial densities that should be positive
    if(ats_only){
      beta.obs <- c(p_y1_1 - p_y1_0)
    }else{
      beta.obs <- c(p_y0_0 - p_y0_1,
                    p_y1_1 - p_y1_0)
      
    }
    return(beta.obs)
  }
  
  
  ## Bootstrap the betas (computed lazily when needed)
  beta.obs_list <- NULL
  get_bootstrap_draws <- function() {
    if (is.null(beta.obs_list)) {
      beta.obs_list <<- compute_bootstrap_draws_clustered(
        f = function(df, d, y, m, ...) {
          get_beta.obs(
            df[[y]],
            df[[d]],
            df[[m]],
            df,
            d,
            yvalues,
            reg_formula = reg_formula
          )
        },
        df = df,
        d = d,
        m = m,
        y = y,
        cluster = cluster,
        numdraws = B,
        fix_n1 = fix_n1,
        return_df = F
      )
    }
    beta.obs_list
  }
  
  ## Get beta.obs using actual data
  beta.obs <- get_beta.obs(yvec, dvec, mvec, df, d, yvalues, reg_formula)
  
  
  ## Get analytic variance for arp and cs
  if (method %in% c("ARP", "CS")) {
    # Get variance matrix of the beta.obs bootsraps
    if (analytic_variance) {
      #Calculate the analytic variance
      sigma.obs <- analytic_variance(yvec = yvec,
                                     dvec = dvec,
                                     mvec = mvec,
                                     df = df,
                                     d = d,
                                     reg_formula = reg_formula,
                                     my_values = my_values,
                                     inequalities_only = TRUE,
                                     clustervec = clustervec,
                                     exploit_binary_m = TRUE)
      
      if (method == "CS" & print_both_var) {
        beta.obs_list <- get_bootstrap_draws()
        sigma.obs_boot <- stats::cov(base::Reduce(base::rbind,
                                                  beta.obs_list))
        print(sigma.obs)
        print(sigma.obs_boot)
      }
    } else {
      beta.obs_list <- get_bootstrap_draws()
      sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                           beta.obs_list))
    }
  }
  
  
  
  ## Run the respective tests
  if (method == "FSST") {
    beta.obs_list <- get_bootstrap_draws()
    # Join beta.obs from actual and boostrapped data
    beta.obs_FSST <- c(list(beta.obs), beta.obs_list)

    # Get variance matrix of the beta.obs bootsraps
    sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                         beta.obs_list))
    
    # Run fsst test
    A.obs <- diag(length(beta.obs))
    A.shp <- matrix(0, ncol = ncol(A.obs))
    A.tgt <- matrix(0, ncol = ncol(A.obs))
    
    beta.shp <- 0
    
    lpm <- lpinfer::lpmodel(A.obs = A.obs,
                            A.shp = A.shp,
                            A.tgt = A.tgt,
                            beta.obs = beta.obs_FSST,
                            beta.shp = beta.shp)
    if (lambda == "dd") {
      lambda <- NA
    } else if (lambda == "ndd") {
      lambda <- 1/sqrt(log(max(length(beta.obs), exp(1))) * log(max(exp(1), log(max(exp(1), n)))))
    }
    
    fsst_result <- lpinfer::fsst(n = n, lpmodel = lpm, beta.tgt = 0, R = B,
                                 weight.matrix = weight.matrix, lambda = lambda)
    
    return(list(result = fsst_result, reject = (fsst_result$pval[1, 2] < alpha)))
  } else if (method == "ARP") {
    if(use_hybrid){
      lf_cv <- HonestDiD:::.compute_least_favorable_cv(X_T = matrix(0,nrow = length(beta.obs)),
                                                       sigma = sigma.obs,
                                                       hybrid_kappa = kappa
      )
      
      hybrid_list <- list(hybrid_kappa = kappa, lf_cv = lf_cv)
      
      arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                 y_T = -beta.obs,
                                                 X_T = matrix(0,nrow = length(beta.obs)),
                                                 sigma = sigma.obs,
                                                 alpha = alpha,
                                                 hybrid_flag = "LF",
                                                 hybrid_list = hybrid_list
      )
    }else{
      arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                 y_T = -beta.obs,
                                                 X_T = matrix(0,nrow = length(beta.obs)),
                                                 sigma = sigma.obs,
                                                 alpha = alpha,
                                                 hybrid_flag = "ARP")
      
    }
    return(arp)
  } else if (method == "CS") {
    #Run Cox and Shi test
    coxandshi <- cox_shi_nonuisance(Y = -beta.obs,
                                    sigma = sigma.obs,
                                    alpha = alpha,
                                    refinement = refinement)
    return(coxandshi)
  }
}
