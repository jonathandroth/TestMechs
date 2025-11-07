#' @title Hypothesis test for the sharp null of full mediation
#' @description This function tests the sharp null that Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. Several options are available for inference.
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable(s)
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param method Method to use. One of "ARP, "CS", "FSST"
#' @param ordering A list with length equal to the cardinality of the support of the mediator variable. The name of each element corresponds to a point in the support, and each element is a vector that collects all m values that are less than or equal to this point. If ordering = NULL, the standard ordering is used. If length(m) > 1, then the default is the elementwise ordering.
#' @param B Bootstrap size, default is 500
#' @param cluster Cluster for bootstrap
#' @param reg_formula Regression formula for the non-experimental setting.
#' @param weight.matrix Weight matrix used to implement FSST. Possible options are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param hybrid_kappa The first-stage size value of the ARP hybrid test. If NULL, the ARP conditional (non-hybrid) test is used. Default is alpha/10
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @param alpha Significance level. Default is 0.05
#' @param rearrange Logical variable indicating whether to rearrange the conditional probabilities to obey monotonicity. Default is FALSE.
#' @param eps_bar Cho and Russell (2023) truncation parameter
#' @param enumerate Enumerate vertices for Cox and Shi (2023) implementation?
#' @param fix_n1 Whether the number of treated units (or clusters) should be fixed in the bootstrap
#' @param lambda Lambda value for FSST. Default is "dd" in which case the "data-driven" lambda recommended by FSST is used. If lambda = "ndd", the non data-driven recommendation is used. See Remark 4.2 of FSST.
#' @param use_nc If the data is clustered, should FSST use the number of clusters for determining lambda (versus total observations). Default is false.
#' @param analytic_variance If TRUE, we use the analytic formula for the variance, rather than a bootstrap. Available if method is ARP or CS. Default is FALSE
#' @param max_defiers_share Bound on the proportion of defiers in the population. Default is 0 which indicates that the monotonicity constraint is imposed.
#' @param new_dof_CS Use the new degrees of freedom formula for Cox and Shi? Default is FALSE.
#' @param use_binary If TRUE, uses ARP and CS implementation that exploits the fact that there are no nuisance parameters when the mediator is binary
#' @param frac_ATs_affected Relaxes the sharp null to test (sum theta_kk nu_k)/sum theta_kk <= frac_ATs_affected. Default is NULL which corresponds to the sharp null.
#' @export
#'
test_sharp_null <- function(df,
                            d,
                            m,
                            y,
                            method,
                            ordering = NULL,
                            B = 500,
                            cluster = NULL,
                            reg_formula = NULL,
                            weight.matrix = "diag",
                            hybrid_kappa = alpha/10,
                            num_Ybins = NULL,
                            alpha = 0.05,
                            rearrange = FALSE,
                            eps_bar = 1e-03,
                            enumerate = FALSE,
                            fix_n1 = TRUE,
                            lambda = "dd",
                            use_nc = FALSE,
                            analytic_variance = ifelse(method == "FSST", FALSE, TRUE),
                            max_defiers_share = 0,
                            new_dof_CS = FALSE,
                            use_binary = NULL,
                            frac_ATs_affected = NULL){

  ## Process the inputted df ----
  # Remove missing values
  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y,
                               reg_formula = if (!is.null(reg_formula)) stats::as.formula(paste(y, reg_formula)) else NULL)


  ## Evaluate whether M is binary
  binary_M <- dplyr::n_distinct(df[m]) == 2

  ## Throw error if M is non-binary but use_binary = TRUE
  if (!is.null(use_binary)) {
    if (use_binary == TRUE & binary_M == FALSE) {
      stop("M is non-binary in data, use_binary cannot be TRUE")
    }
  }
  else if (is.null(use_binary)) {
    use_binary = binary_M
  }

  ## Set frac_ATs_affected to NULL for now, will allow for other values later on
  ## frac_ATs_affected <- NULL

  ## Do not run the binary test when frac_ATs_affected is not null, since we want to test for the fraction of always-takers
  if (!is.null(frac_ATs_affected) | !(max_defiers_share == 0)) {
    if (method == "toru") {
      stop("You cannot allow for defiers when method == `toru`")
    }
    use_binary = FALSE
  }

  ## Make sure lpinfer is installed if method == "FSST"
  if (method == "FSST" & !require(lpinfer)) {
    stop("You must install package lpinfer to use the FSST method. You can install it by running
devtools::install_github('conroylau/lpinfer')")
  }

  ## Print confirmation, can delete this if you prefer
  ## message(paste0("It is ", binary_M, " that M in data is binary"))
  ## message(paste0("It is ", use_binary, " that binary test is used"))

  ## Pass to more efficient algorithm when M is binary
  if (use_binary) {
    if (method %in% c("ARP", "CS", "FSST")) {
      result <- test_sharp_null_binary_m(df,
                                         d,
                                         m,
                                         y,
                                         method = method,
                                         reg_formula = reg_formula,
                                         ordering = ordering,
                                         B = B,
                                         cluster = cluster,
                                         weight.matrix = weight.matrix,
                                         alpha = alpha,
                                         kappa = hybrid_kappa,
                                         use_hybrid = T,
                                         num_Ybins = num_Ybins,
                                         fix_n1 = fix_n1,
                                         lambda = lambda,
                                         analytic_variance = analytic_variance,
                                         refinement = TRUE)
      return(result)
    }
    else if (method == "toru") {
      if(!is.null(reg_formula)){
        stop("reg_formula should be NULL if method = toru")
      }else{
      result <- test_sharp_null_toru(df, d, m, y, B = B, alpha = alpha,
                                     num_Ybins = NULL, cluster = cluster)

      return(result)
      }
    }
    else {
      stop("Method must be either ARP, CS, FSST or toru if use_binary = TRUE")
    }
  }

  # Sample size
  n <- nrow(df)

  # If M is multivariate convert to univariate and use elementwise ordering
  if (length(m) > 1) {

    m_supp <- unique(df[m])
    ordering <- vector("list", length = nrow(m_supp))
    names(ordering) <- 1:nrow(m_supp)

    # Create ordering and replace multivariate M by a univariate M
    uni_m <- numeric(n)
    for (k in 1:nrow(m_supp)) {
      m_k <- m_supp[k,]
      ordering[[as.character(k)]] <-
        which(apply(m_supp, 1, function(m_supp_elem) all(m_supp_elem <= m_k)))
      uni_m[colSums(as.vector(m_k) == t(df[m])) == 2] <- k
    }

    df$m <- uni_m
    m <- "m"
  }

  #Discretize y if needed
  if(!is.null(num_Ybins)){
    df[[y]] <- discretize_y(yvec = df[[y]], numBins = num_Ybins)
  } else {
    continuous_y_flag <- n / length(unique(df[[y]])) <= 30
    if (continuous_y_flag) {
      message("Y variable might be continuous. Discretize it by specifying num_Ybins. Default num_Ybins = 5 is used now.")
      num_Ybins <- 5
      df[[y]] <- discretize_y(yvec = df[[y]], numBins = num_Ybins)
    }
  }

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  if(is.null(cluster)){
    clustervec <- 1:length(yvec)
  }else{
    clustervec <- df[[cluster]]
  }

  ## Construct the A matrices and beta.shp ----

  #Specify whether method requires us to input only inequalities
  inequalities_only <- ifelse(method %in% c("ARP","CS"),
                              TRUE, FALSE )

  #Construct the relevant A matrices and beta.shp
  A_list <- construct_Aobs_Ashp_betashp(yvec = yvec,
                                        mvec = mvec,
                                        ordering = ordering,
                                        inequalities_only = inequalities_only,
                                        max_defiers_share = max_defiers_share,
                                        frac_ATs_affected = frac_ATs_affected)

  A.shp <- A_list$A.shp
  A.obs <- A_list$A.obs
  beta.shp <- A_list$beta.shp

  ## Construct beta.obs and beta.obs_list ----

  #The following values are used by beta.obs
  #Note we need to do this here so that the yvalues remain constant across boostrap draws
  yvalues <- unique(yvec)
  mvalues <- unique(mvec)
  my_values <- purrr::cross_df(list(m=mvalues,y=yvalues)) %>%
    dplyr::arrange(m,y) %>%
    dplyr::select(y,m)

  # Override analytic var if method = FSST
  if (method == "FSST") {
    analytic_variance <- FALSE
  }

  # Compute beta.obs at the observed data
  beta.obs <- get_beta.obs_fn(yvec = yvec,
                              dvec = dvec,
                              mvec = mvec,
                              df = df,
                              d = d,
                              reg_formula = reg_formula,
                              inequalities_only = inequalities_only,
                              yvalues = yvalues,
                              mvalues = mvalues,
                              my_values = my_values,
                              rearrange = rearrange)

  if(!analytic_variance | method %in% c("FSST")){
    #Bootstrap to get beta.obs_list
    beta.obs_list <- compute_bootstrap_draws_clustered(
      f =
        function(df,d,y,m,...){get_beta.obs_fn(
          yvec = df[[y]],
          dvec = df[[d]],
          mvec = df[[m]],
          df = df,
          d = d,
          reg_formula = reg_formula,
          inequalities_only = inequalities_only,
          yvalues = yvalues,
          mvalues = mvalues,
          my_values = my_values,
          rearrange = rearrange)},
      df = df,
      d = d,
      m = m,
      y = y,
      cluster = cluster,
      numdraws = B,
      return_df = F,
      fix_n1 = fix_n1)
  }

  if (analytic_variance){
    #Calculate the analytic variance
    sigma.obs <- analytic_variance(yvec = yvec,
                                   dvec = dvec,
                                   mvec = mvec,
                                   df = df,
                                   d = d,
                                   reg_formula = reg_formula,
                                   my_values = my_values,
                                   inequalities_only = inequalities_only,
                                   clustervec = clustervec)
  }
  ## Pass to the relevant moment inequality procedure ----

  if(method %in% c("FSST")){
    # Define target parameter
    A.tgt <- A_list$A.tgt

    # Run FSST
    if (analytic_variance) {
      beta.obs_FSST <- list(c(list(beta.obs),
                              beta.obs_list),
                            sigma.obs)
    } else {
      beta.obs_FSST <- c(list(beta.obs), beta.obs_list)
    }

    lpm <- lpinfer::lpmodel(A.obs = A.obs,
                            A.shp = A.shp,
                            A.tgt = A.tgt,
                            beta.obs = beta.obs_FSST,
                            beta.shp = beta.shp)

    if(is.null(cluster) | use_nc == FALSE){
      n <- NROW(df)
    }else{
      n <- length(unique(df[[cluster]]))
    }

    if (lambda == "dd") {
      lambda <- NA
    } else if (lambda == "ndd") {
      lambda <- 1/sqrt(log(max(length(beta.obs), exp(1))) * log(max(exp(1), log(max(exp(1), n)))))
    }

    fsst_result <- lpinfer::fsst(n = n, lpmodel = lpm, beta.tgt = 0, R = B,
                                 weight.matrix = weight.matrix, lambda = lambda)

    return(list(result = fsst_result, reject = (fsst_result$pval[1, 2] < alpha)))
  }


  if(method %in% c("ARP", "CS")){

    if(!analytic_variance){
      # Get variance matrix of the beta.obs boostraps
      sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                           beta.obs_list))
    }

    #Add the shape constraints as moments
    beta <- c(beta.obs, beta.shp)
    A <- rbind(A.obs, A.shp)

    #Update the covariance matrix to have zero in blocks for the shape constraints
    sigma <- rbind(cbind(sigma.obs, matrix(0,
                                           nrow = NROW(sigma.obs),
                                           ncol = NROW(A.shp) ) ),
                   matrix(0, nrow = NROW(A.shp), ncol = NCOL(sigma.obs) + NROW(A.shp) ))

    #Run the relevant test
    if (method == "ARP") {

      #Test if sigma is numerically not psd. If so, add small amount of noise
      min_eig <- base::min(base::eigen(sigma, only.values = TRUE)$values)
      if (min_eig < 0) {
        sigma <- sigma + diag(10 * abs(min_eig),
                              ncol = NROW(sigma),
                              nrow = NROW(sigma))
      }

      if (is.null(hybrid_kappa)) {
        arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                   y_T = beta,
                                                   X_T = A,
                                                   sigma = sigma,
                                                   alpha = alpha,
                                                   hybrid_flag = "ARP")
      } else {
        lf_cv <- HonestDiD:::.compute_least_favorable_cv(X_T = A,
                                                         sigma = sigma,
                                                         hybrid_kappa = hybrid_kappa)

        hybrid_list <- list(hybrid_kappa = hybrid_kappa,
                            lf_cv = lf_cv)

        arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                   y_T = beta,
                                                   X_T = A,
                                                   sigma = sigma,
                                                   alpha = alpha,
                                                   hybrid_flag = "LF",
                                                   hybrid_list = hybrid_list)
      }

      return(arp)
    }

    if(method == "CS"){

      # We have beta - C_Z delta <= d_Z
      # C_Z = A; d_Z = 0
      # We want to define B_Z and m so that var(m) has full rank
      # This can be done by taking the rows corresponing to postive eigenvals

      ## Later, convert in C * delta <= m for to make use of the replication package
      C_Z <- A
      d_Z <- 0


      if (min(base::eigen(sigma, only.values = T)$values) < 1e-08){

        #If sigma is not full-rank we extract the full rank component
        eigenvecs <- eigen(sigma)$vectors
        eigenvals <- eigen(sigma)$values

        tol <- 1e-08
        positive_indices <- which(eigenvals > tol)
        ## qr(sigma)$rank

        B_Z <- eigenvecs[,positive_indices, drop = FALSE]

        # Get  "reduced" beta
        # By constuction, B_Z %*% t(B_Z) %*% beta - beta = constant
        beta_red <- t(B_Z) %*% beta

        V_red <- diag(eigenvals[positive_indices]) # variance of reduced beta

        # Adjust the linear programming parameters accordingly
        # We have B_Z %*% beta.red - C_Z delta <= d_Z + B_Z %*% beta.red - beta, following CS notation
        d_Z <- d_Z + B_Z %*% beta_red - beta
        sigma <- V_red
        sigmaInv <- solve(sigma)
      } else {
        beta_red <- beta.obs
        B_Z_red <- diag(lengt(beta.obs))
        sigmaInv <- solve(sigma)
      }

      d_nuis <- ncol(A)  # number of nuisance parameters
      d_ineq <- nrow(A)  # number of inequalities

      # Perform CC Test

      # Finding T (Test Statistic)
      # Amat x <= d_Z
      Amat <- rbind(cbind(B_Z, - C_Z),
                    cbind(matrix(0, nrow = d_nuis, ncol = ncol(B_Z)), diag(d_nuis)))
      u <- c(d_Z, rep(1, d_nuis))
      l <- c(rep(-Inf, length(d_Z)), rep(0, d_nuis))
      Dmat <- matrix(0, nrow = (nrow(sigma) + d_nuis), ncol = (nrow(sigma) + d_nuis))
      Dmat[1:nrow(sigma), 1:nrow(sigma)] <- 2 * sigmaInv
      dvec <- matrix(c(2 * sigmaInv %*% beta_red, numeric(d_nuis)), ncol = 1)


      osqp_tol <- 1e-8 # important - smaller tol level don't work well
      osqp_opts <- list(verbose = FALSE,
                        eps_abs = osqp_tol,
                        eps_rel = osqp_tol)

      qp <- osqp::solve_osqp(P = Dmat, q = -dvec, A = Amat, l = l, u = u, pars = osqp_opts)
      ## qp <- quadprog::solve.QP(Dmat = Dmat,
      ##                          dvec = dvec,
      ##                          Amat = Amat,
      ##                          bvec = - d_Z)



      ## beta.obs_red_star <- qp$solution[1:nrow(sigma)]
      beta_red_star <- qp$x[1:nrow(sigma)]
      T_CC <- t(beta_red - beta_red_star) %*%
        sigmaInv %*%
        (beta_red - beta_red_star)

      ## Avoids some numerical instability issues
      beta_red_star[abs(beta_red_star) < osqp_tol] <- 0

      if (new_dof_CS) {

        Amat_aug <- rbind(Amat,
                          -cbind(matrix(0, nrow = d_nuis, ncol = ncol(B_Z)), diag(d_nuis)))
        u_aug <- c(u, rep(0, d_nuis))

        Khat <- which(abs(u_aug - Amat_aug %*% qp$x) < osqp_tol^(1/2))

        ## I_Khat <- diag(ncol(Amat_aug))[Khat, , drop = FALSE]

        qr_tol <- 1e-5 # important; don't be too accurate

        dof_n  <- qr(Amat_aug[Khat, ], tol = qr_tol)$rank -
                             qr(Amat_aug[Khat, -(1:ncol(B_Z))], tol = qr_tol)$rank
        cv_CC <- qchisq(1 - alpha, df = dof_n)

        return(list(T_CC = T_CC,
                    cv_CC = cv_CC,
                    df = dof_n,
                    reject = (T_CC > cv_CC),
                    pval = 1-stats::pchisq(q = T_CC, df = dof_n)))
      }

      # Equivalent expression
      ## T_CC <- qp$info$obj_val + t(beta_red) %*% sigmaInv %*% beta_red

      # Find the Degree of Freedom:
      ## tol <- 1e-6

      ## if (enumerate) {
      ##   A_vert <- rbind(-diag(d_ineq),
      ##                   t(C_Z),
      ##                   -t(C_Z),
      ##                   rep(1, d_ineq),
      ##                   -rep(1, d_ineq))
      ##   b_vert <- c(numeric(d_ineq), 2 * numeric(d_nuis), 1, -1)
      ##   H <- vertexenum::enumerate.vertices(A = A_vert, b = b_vert)
      ##   A_Z <- H %*% B_Z_red
      ##   b_Z <- H %*% d_Z

      ##   # Binding constraints

      ##   dof_n <- sum(abs(A_Z %*% beta.obs_red_star - b_Z) > tol)
      ## } else {
      # Define A_ineq0 and b_ineq0
      ## A_ineq0 <- - diag(d_ineq)
      I_d_ineq <- diag(d_ineq)
      ## b_ineq0 <- matrix(0, nrow = d_ineq, ncol = 1)

      ## Finding I_J0 as in Section A.3.1 of Cox and Shi (2023)

      model <- list()
      model$A <- rbind(t(C_Z),  rep(1, d_ineq))
      model$obj <- d_Z - B_Z %*% beta_red_star
      model$modelsense <- "min"
      model$rhs <- c(rep(0, d_nuis), 1)
      model$sense <- rep("=", nrow(model$A))
      model$lb <- rep(0, d_ineq)

      ## gurobi_tol <- 1e-7
      ## params <- list(OutputFlag=0, OptimalityTol = gurobi_tol, FeasibilityTol = gurobi_tol)
      ## result <- gurobi::gurobi(model, params)

      result <- Rglpk::Rglpk_solve_LP(model$obj, model$A, rep("==", nrow(model$A)), model$rhs)
      result$status <- ifelse(result$status == 0, "OPTIMAL", "ELSE")
      result$objval <- result$optimum

      ## result$x

      V_min <- result$objval
      flag <- result$status

      if ((flag != "OPTIMAL") | (V_min >= 0.00005)) {

        dof_n <- 0
        cv_CC <- qchisq(1 - alpha, df = dof_n)

        return(list(T_CC = T_CC, cv_CC = cv_CC, df = dof_n, reject = (T_CC > cv_CC)))
      }

      psis <- rep(NA, d_ineq)

      for (j in 1:d_ineq) {

        model <- list()
        model$A <- rbind(t(C_Z), t(d_Z - B_Z %*% beta_red_star), rep(1, d_ineq))
        model$obj <- -I_d_ineq[j, ]
        model$modelsense <- "min"
        model$rhs <- c(rep(0, d_nuis), V_min, 1)
        ## model$rhs <- c(rep(0, d_nuis), 0, 1)
        ## model$sense <- rep("=", nrow(model$A))
        model$sense <- rep("=", nrow(model$A))
        model$lb <- rep(0, d_ineq)

        result <- Rglpk::Rglpk_solve_LP(model$obj, model$A, rep("==", nrow(model$A)), model$rhs)
        result$status <- ifelse(result$status == 0, "OPTIMAL", "ELSE")
        result$objval <- result$optimum

        ## params <- list(OutputFlag=0)
        ## result <- gurobi::gurobi(model, params)

        if (result$status != "OPTIMAL") {
          psis[j] <- -1
        } else {
          psis[j] <- result$objval
        }

      }

      # Define A_eq0, A_eq, and b_eq
      # mstar = - d_Z + B_Z %*% beta_red*
      ## mstar <- d_Z - B_Z %*% beta_red_star
      ## A_eq0 <- rbind(-t(C_Z), t(mstar), rep(1, d_ineq))
      ## A_eq <- rbind(-t(C_Z), rep(1, d_ineq))
      ## b_eq <- c(rep(0, d_nuis), 1)

      # Perform linear programming to find Vmu_min
      ## lp_options <- lpSolve::lp.control(display = "silent", simplex = "dual")
      ## result <- lpSolve::lp("min", mstar, rbind(A_ineq0, A_eq),
      ##                       c(rep(">=", d_ineq), rep("=", nrow(A_eq))),
      ##                       c(b_ineq0, b_eq))
      ## Vmu_min <- result$objval
      ## flag <- result$status

      ## model <- list()
      ## model$A <- A_eq
      ## model$obj <- mstar
      ## model$modelsense <- "min"
      ## model$rhs <- b_eq
      ## model$sense <- rep("=", nrow(A_eq))
      ## model$lb <- rep(0, d_ineq)

      ## params <- list(OutputFlag=0)
      ## result <- gurobi::gurobi(model, params)

      ## Vmu_min <- result$objval
      ## flag <- result$status

      # Check conditions and calculate dof_n
      ## if (flag != "OPTIMAL") {
      ##   dof_n <- 0
      ##   ## } else if (flag == 2) { ## for lpSolve
      ## } else if (Vmu_min >= 0.00005) {
      ##   dof_n <- 0
      ## } else {
      ##   # Initialize nb_ineq_min
      ##   nb_ineq_min <- rep(NA, d_ineq)
      ##   counter <- 1

      ##   # Iterate through each j
      ##   for (bj in 1:d_ineq) {
      ##     # Find the largest b_j allowed

      ##     ## result_j <- lpSolve::lp("min", A_ineq0[bj,],
      ##     ##                         rbind(A_eq0, diag(d_ineq), diag(d_ineq)),
      ##     ##                         c(rep("=", nrow(A_eq0)), rep(">=", d_ineq), rep("<=", d_ineq)),
      ##     ##                         c(rep(0, d_nuis), Vmu_min, 1, b_ineq0, b_ineq0 + 1))

      ##     model <- list()
      ##     model$A <- A_eq0
      ##     model$obj <- A_ineq0[bj, ]
      ##     model$modelsense <- "min"
      ##     model$rhs <- c(rep(0, d_nuis), Vmu_min, 1)
      ##     model$sense <- rep("=", nrow(A_eq0))
      ##     model$lb <-  0
      ##     model$ub <- 1

      ##     result_j <- gurobi::gurobi(model, params)
      ##     # Update nb_ineq_min
      ##     ## if (result_j$status == 2) {# for lpSolve
      ##     if (result_j$status != "OPTIMAL") {
      ##       nb_ineq_min[bj] <- - 1
      ##       counter <- counter + 1
      ##     } else {
      ##       nb_ineq_min[bj] <- result_j$objval
      ##     }
      ##   }

      # Collect rows of A_ineq corresponding to implicit equalities
      ## A_impeq <- A_ineq0[(0 - nb_ineq_min) < tol, ]

      psi_tol <- 1e-6

      A_impeq <- I_d_ineq[ (-psis) < psi_tol, , drop = FALSE]

      A_full_eq <- rbind(A_impeq, -t(C_Z), t(B_Z %*% beta_red_star - d_Z))

      # Calculate rank of A_full_eq
      qr_tol <- 1e-5 # important; don't be too accurate
      rkA <- qr(t(A_full_eq) %*% A_full_eq, tol = qr_tol)$rank

      if (qr(B_Z)$rank == d_ineq) {

        # Combine A_impeq and A_eq0
        ## A_full_eq <- rbind(A_impeq, rbind(-t(C_Z), t(B_Z %*% beta_red_star - d_Z)))

        ## # Calculate rank of A_full_eq
        ## rkA <- qr(A_full_eq)$rank

        # Calculate dof_n
        dof_n <- d_ineq - rkA

      } else if (qr(B_Z)$rank < d_ineq) {

        ## G <- cbind(t(A_impeq), C_Z, B_Z %*% beta_red_star - d_Z)
        ## qr_tG <- qr(t(G))
        ## rank_G <- qr_tG$rank

        G <- t(A_full_eq)
        rank_G <- rkA

        if (rank_G == nrow(G)) {

          dof_n <- 0

        } else {

          ## G1_ind <- qr_tG$pivot[1:rank_G]
          ## G2_ind <- setdiff(1:nrow(G), G1_ind)

          ## G1 <- G[G1_ind, , drop = FALSE]
          ## G2 <- G[G2_ind, , drop = FALSE]

          ## B_Z1 <- B_Z[G1_ind, , drop = FALSE]
          ## B_Z2 <- B_Z[G2_ind, , drop = FALSE]

          ## Gamma <- - solve(G1 %*% t(G1)) %*% G1 %*% t(G2)
          ## dof_n <- qr(t(Gamma) %*% B_Z1 + B_Z2)$rank

          G1_ind <- qr(t(G), tol = qr_tol)$pivot[1:rank_G]
          G2_ind <- setdiff(1:nrow(G), G1_ind)

          G1 <- G[G1_ind, , drop = FALSE]
          G2 <- G[G2_ind, , drop = FALSE]

          B_Z1 <- B_Z[G1_ind, , drop = FALSE]
          B_Z2 <- B_Z[G2_ind, , drop = FALSE]

          ## Gamma <- - solve(G1 %*% t(G1), tol = .Machine$double.eps^3) %*% G1 %*% t(G2)
          Gamma <- - qr.coef(qr(t(G1), tol = qr_tol), t(G2))
          dof_n <- qr(t(Gamma) %*% B_Z1 + B_Z2, tol = qr_tol)$rank

        }
      }

      ## }
      # Enumerating the vertices using (outdated) package vertexenum
      # Finds vertices for the polytope Ax <= b

      cv_CC <- qchisq(1 - alpha, df = dof_n)
      return(list(T_CC = T_CC,
                  cv_CC = cv_CC,
                  df = dof_n,
                  reject = (T_CC > cv_CC),
                  pval = 1-stats::pchisq(q = T_CC, df = dof_n)))
    }
  }


  stop("method must be one of ARP, CS, FSST")

}


#Function for creating the matrices A.obs and A.shp and the vector beta.shp
construct_Aobs_Ashp_betashp <- function(yvec,
                                        mvec,
                                        ordering,
                                        inequalities_only = F,
                                        max_defiers_share = 0,
                                        frac_ATs_affected = NULL){

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

  # If testing the fraction of ATs, then:
  # x = (theta, delta, zeta, kappa, eta, iota, gamma, epsilon)
  # iota: defined to be v_k theta_kk
  # gamma: nuisance par to convert inequalities to equalities
  # epsilon: nuisance par to create target par

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "zeta" = K,
                   "kappa" = d_y * K, "eta" = K)

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

  # If testing the fraction of ATs, include the additional params
  if (!is.null(frac_ATs_affected)) {

    # iota is eta in Jon's note, for implementing the fraction of always takers affected.
    # gamma is nuisance par to convert inequalities to equalities.
    # epsilon: nuisance par to create target par for the new constraints

    par_lengths <- c(par_lengths, "iota" = K, "gamma" = K + 1, "epsilon" = K + 1)
    parameter_types <- c(parameter_types,
                         rep("iota", par_lengths[["iota"]]),
                         rep("gamma", par_lengths[["gamma"]]),
                         rep("epsilon", par_lengths[["epsilon"]]))
    iota_indices <- which(parameter_types == "iota")
    gamma_indices <- which(parameter_types == "gamma")
    epsilon_indices <- which(parameter_types == "epsilon")
  }

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
  if (max_defiers_share == 0) {
    A.shp <- matrix(0, nrow = K, ncol = len_x)
  } else {
    A.shp <- matrix(0, nrow = K + 1, ncol = len_x)
  }


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


  ## If testing the fraction of ATs, need to update previously defined inequalities and include new inequalities
  ## For more details, please refer to Jon's notes.
  if (!is.null(frac_ATs_affected)) {

    ## Coefficients for iota, for equation 2
    A.shp[1:K, sum(par_lengths[1:5]) + 1:K] <- diag(K)

    ## For equation 3 and 4
    temp <- matrix(0, nrow = K + 1, ncol = len_x)

    # for equation 4
    for (k in 1:K) {
      temp[k, (k-1) * K + k] <- 1    # theta_kk has coefficient 1
      temp[k, sum(par_lengths[1:5]) + k] <- -1  # iota_k has coefficient -1
      temp[k, sum(par_lengths[1:6]) + k] <- -1  # gamma is slack parameter for turning ineq into eq
      temp[k, sum(par_lengths[1:7]) + k] <- 1 # epsilon is the FSST target param
    }

    # for equation 3
    temp[K+1, 1:par_lengths[1]] <- colSums(temp[1:K, 1:par_lengths[1]]) * frac_ATs_affected  # summing thetas
    temp[K+1, sum(par_lengths[1:5]) + 1:K] <- -1  # summing iotas
    temp[K+1, sum(par_lengths[1:7])] <- -1 # gamma converting inequalities into equalities
    temp[K+1, sum(par_lengths[1:7]) + K + 1] <- 1 # epsilon is the FSST target param

    A.shp <- rbind(A.shp, temp)
  }


  if (max_defiers_share == 0) {
    if (inequalities_only) {
      # Because these are for the inequalities to become equalities. eta is the FSST target param which is a nuisance param
      if (is.null(frac_ATs_affected)) {
        A.shp <- A.shp[, -c(l_gt_k_inds, zeta_indices, kappa_indices, eta_indices)]
      } else { # If testing the fraction of ATs, remove the additional nuisance pars
        A.shp <- A.shp[, -c(l_gt_k_inds, zeta_indices, kappa_indices, eta_indices, gamma_indices, epsilon_indices)]
      }

      #Add shape constraint that all parameters are >= 0 (this is not enforced by ARP)
      A.shp <- rbind(A.shp, diag(NCOL(A.shp)))

    } else {
      #Remove the extraneous thetas only
      A.shp <- A.shp[, -l_gt_k_inds]
    }

  } else {

    A.shp[K + 1, l_gt_k_inds] <- - 1

    if (inequalities_only) {
      #Remove both the extraneous thetas *and* kappa,zeta,eta
      if (is.null(frac_ATs_affected)) {
        A.shp <- A.shp[, -c(zeta_indices, kappa_indices, eta_indices)]
      } else { # If testing the fraction of ATs, remove the additional nuisance pars
        A.shp <- A.shp[, -c(zeta_indices, kappa_indices, eta_indices, gamma_indices, epsilon_indices)]
      }
      #Add shape constraint that all parameters are >= 0 (this is not enforced by ARP)
      A.shp <- rbind(A.shp, diag(NCOL(A.shp)))

    } else {

      # Inequalities to equalities for max_defiers_share
      if (!is.null(frac_ATs_affected)) {
        A.shp <- cbind(A.shp, c( rep(0, K), -1, rep(0, K+1) )) # A slackness param which is greater than or equal to 0 to impose equality
      } else {
        A.shp <- cbind(A.shp, c(rep(0, K), -1)) # A slackness param which is greater than or equal to 0 to impose equality
      }

      # Update the list of parameters
      par_lengths <- c(par_lengths, "defier_zeta" = 1)
      parameter_types <- c(parameter_types, "defier_zeta")
      defier_zeta_indices <- which(parameter_types == "defier_zeta")

      len_x <- sum(par_lengths)
    }
  }

  beta.shp <- rep(0, NROW(A.shp))

  if (max_defiers_share != 0) {
    beta.shp[K + 1]  <- - max_defiers_share
  }


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
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- -1
  }

  if (max_defiers_share == 0) {
    if (inequalities_only) {
      # Remove both the extraneous thetas *and* kappa,zeta,eta
      if (is.null(frac_ATs_affected)) {
        A.obs <- A.obs[, -c(l_gt_k_inds, zeta_indices, kappa_indices, eta_indices)]
      } else { # If testing the fraction of ATs, remove the additional nuisance pars
        A.obs <- A.obs[, -c(l_gt_k_inds, zeta_indices, kappa_indices, eta_indices, gamma_indices, epsilon_indices)]
      }

      # The first 2K rows of A.obs are equality constraints, so we duplicate
      # them with opposite signs to get equalities as inequalities
      A.obs <- rbind(A.obs[1:(2 * K),],
                     -A.obs[1:(2 * K),],
                     A.obs[(2 * K + 1):NROW(A.obs), ])
    } else {
      #Remove the extraneous thetas only
      A.obs <- A.obs[, -c(l_gt_k_inds)]
    }

  } else {
    if (inequalities_only) {
      # Remove kappa,zeta,eta
      if (is.null(frac_ATs_affected)) {
        A.obs <- A.obs[, -c(zeta_indices, kappa_indices, eta_indices)]
      } else { # If testing the fraction of ATs, remove the additional nuisance pars
        A.obs <- A.obs[, -c(zeta_indices, kappa_indices, eta_indices, gamma_indices, epsilon_indices)]
      }

      # The first 2K rows of A.obs are equality constraints, so we duplicate
      # them with opposite signs to get equalities as inequalities
      A.obs <- rbind(A.obs[1:(2 * K),],
                     -A.obs[1:(2 * K),],
                     A.obs[(2 * K + 1):NROW(A.obs), ])
    }
  }

  # Create A.tgt (only used for lpinfer functions)
  A.tgt <- numeric(len_x)
  A.tgt[sum(par_lengths[1:4]) + (1:K)] <- 1

  # If testing the fraction of ATs, include the additional target pars
  if (!is.null(frac_ATs_affected)) {
    A.tgt[sum(par_lengths[1:7]) + 1:(K+1)] <- 1
  }

  # Remove the defier theta column
  if (max_defiers_share == 0) {
    A.tgt <- A.tgt[-l_gt_k_inds]
  }



  return(list(A.obs = A.obs,
              A.shp = A.shp,
              A.tgt = A.tgt,
              beta.shp = beta.shp,
              par_length = length(A.tgt)))
}



# Constructing beta_obs
get_beta.obs_fn <- function(yvec, dvec, mvec, df, d, reg_formula = NULL, inequalities_only,
                            yvalues, mvalues, my_values, rearrange = FALSE) {
  #Get frequencies for all possible values of (y,m) | D=0
  if(!is.null(reg_formula)){
    #Use regression approach
    #Parse reg_formula to check for IV or covariates
    #iv_spec <- extract_iv(reg_formula, d)

    #XX should do some basic sanity checks on the regression spec here

    p_ym_0_vec <- numeric(nrow(my_values))

      for (j in 1:nrow(my_values)) {

        # Create indicator for current (y,m) combo
        df$lhs <- as.numeric(yvec == my_values$y[j] & mvec == my_values$m[j])

        #Interact indicator with (D-1)
        # This is because treatment effect on (D-1)*Y is simply Y(0)
        df$lhs <- (df[[d]]-1) * df$lhs

        #Check if LHS is constant (since if so, fixest throws an error)
        #If so, we difference btwn treatmetn and control is zero, so report 0
        if(var(df$lhs) == 0){
          p_ym_0_vec[j] <- 0
          next
        }

        fml <- stats::as.formula( base::paste("lhs", reg_formula) )

        reg <- fixest::feols(fml = fml, data = df)

        if( d %in% rownames(reg$coeftable) ){
          p_ym_0_vec[j] <- reg$coeftable[d,"Estimate"]
        }else if( base::max( base::grepl(pattern = base::paste0(d," ="), x = rownames(reg$coeftable)) ) >0 ){
          p_ym_0_vec[j] <- reg$coeftable[which(base::grepl(base::paste0(d," ="),x = rownames(reg$coeftable) )), "Estimate"]
        }else{
            stop("The treatment variable does not appear to be on the RHS of the provided reg_formula")
          }

    }
  } else {
      # Use original frequency approach (randomized D)
      p_ym_0_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                                   .f = ~mean(yvec[dvec == 0] == my_values$y[.x]
                                              & mvec[dvec == 0] == my_values$m[.x]) )

    }
  #Rearrange into a matrix where rows correspond to values of y and cols vals of m
  p_ym_d0 <-
    cbind(my_values, p_ym_0_vec) %>%
    reshape2::dcast(y ~ m, value.var = "p_ym_0_vec") %>%
    dplyr::select(-y) %>%
    as.matrix()


  #Get frequencies for all possible values of (y,m) | D=1
   if(!is.null(reg_formula)){
    #Use regression approach
    #Parse reg_formula to check for IV or covariates
    #iv_spec <- extract_iv(reg_formula, d)
    p_ym_1_vec <- numeric(nrow(my_values))

    for (j in 1:nrow(my_values)) {
      # Create indicator for current (y,m) combo
      df$lhs <- as.numeric(yvec == my_values$y[j] & mvec == my_values$m[j])

      #Interact indicator with D
      # This is because treatment effect on D*Y is simply Y(1)
      df$lhs <- df[[d]] * df$lhs

      #Check if LHS is constant (since if so, fixest throws an error)
      #If so, we difference btwn treatmetn and control is zero, so report 0
      if(var(df$lhs) == 0){
        p_ym_1_vec[j] <- 0
        next
      }

      fml <- stats::as.formula( base::paste("lhs", reg_formula) )

      reg <- fixest::feols(fml = fml, data = df)

      if( d %in% rownames(reg$coeftable) ){
        p_ym_1_vec[j] <- reg$coeftable[d,"Estimate"]
      }else if( base::max( base::grepl(pattern = base::paste0(d," ="), x = rownames(reg$coeftable)) ) >0 ){
        p_ym_1_vec[j] <- reg$coeftable[which(base::grepl(base::paste0(d," ="),x = rownames(reg$coeftable) )), "Estimate"]
      }else{
        stop("The treatment variable does not appear to be on the RHS of the provided reg_formula")
      }

      }} else{
    # Use original frequency approach (randomized D)
    p_ym_1_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                                 .f = ~mean(yvec[dvec == 1] == my_values$y[.x]
                                            & mvec[dvec == 1] == my_values$m[.x]) )
 }
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

  if (rearrange) {

    cdf_m_d0 <- pmax(cumsum(p_m_d0), cumsum(p_m_d1))
    cdf_m_d1 <- pmin(cumsum(p_m_d0), cumsum(p_m_d1))

    p_m_d0 <- diff(c(0, cdf_m_d0))
    p_m_d1 <- diff(c(0, cdf_m_d1))

  }

  if (inequalities_only) {
    #Duplicate the first two sets of rows with opposite signs
    # to cast equalities as inequalities
    beta.obs <- c(p_m_d0, p_m_d1,
                  -p_m_d0, -p_m_d1,
                  p_ym_1_vec - p_ym_0_vec)
  } else {
    beta.obs <- c(p_m_d0, p_m_d1,
                  p_ym_1_vec - p_ym_0_vec)
  }

  return(beta.obs)
}

# parse OLS/IV regression formula for fixest
extract_iv <- function(reg_formula, d){
  reg_str <- if (inherits(reg_formula, "formula")) {
    paste(deparse(reg_formula), collapse = " ")
  } else {
    as.character(reg_formula)
  }
  reg_str <- trimws(sub("^~", "", reg_str))  # drop leading "~"

  # split by pipe if user left FE / etc (we ignore tail)
  split_pipe <- strsplit(reg_str, "\\|", fixed = FALSE)[[1]]
  rhs_main   <- trimws(split_pipe[1])

  # Initialize output
  out <- list(is_iv    = FALSE,
              treat    = d,
              instr    = character(0),
              controls = character(0),
              added_treat  = FALSE)

  # detect "( ... = ... )"
  if (grepl("\\([^)]*=[^)]*\\)", rhs_main)) {
    # IV branch
    out$is_iv <- TRUE
    iv_part <- sub(".*\\(([^)]*)\\).*", "\\1", rhs_main)
    sides   <- strsplit(iv_part, "=", fixed = TRUE)[[1]]

    out$treat <- trimws(unlist(strsplit(sides[1], "+", fixed = TRUE)))
    out$instr <- trimws(unlist(strsplit(sides[2], "+", fixed = TRUE)))

    rhs_controls <- gsub("\\([^)]*\\)", "", rhs_main)
    ctrls_raw    <- trimws(unlist(strsplit(rhs_controls, "+", fixed = TRUE)))
    out$controls <- setdiff(ctrls_raw, c(out$treat, ""))

    # Enforce that the IV endogenous variable equals d
    if (!all(out$treat == d)){
      stop(
        "In IV syntax, the endogenous treatment inside '(...=...)' must equal d = '",
        d, "'."
      )
    }

  } else {
    # OLS branch
    vars <- trimws(unlist(strsplit(rhs_main, "+", fixed = TRUE)))
    vars <- vars[nzchar(vars)]
    out$controls <- setdiff(vars, d)

    # If user DID NOT include the treatment, we will auto-add it AND warn
    if (!any(vars == d)){
      out$added_treat <- TRUE
      warning(
        "The treatment variable '", d, "' was not found in the provided reg_formula;",
        "I have added it as a regressor. Please edit reg_formula if that was not your intention."
      )
    }
  }
  # Final sanity check
  if (!d %in% c(out$treat, out$controls)) {
    stop("Treatment variable '", d, "' not found in reg_formula. Please include it.")
  }
  out
}


#Function to get the IFs for beta_obs and its subcomponents
get_IFs <- function(yvec, dvec, mvec, df, d, reg_formula = NULL, my_values, mvalues = unique(my_values$m),
                    inequalities_only = T, exploit_binary_m = FALSE){

  n <- length(yvec)
  n0 <- sum(dvec == 0)
  n1 <- sum(dvec == 1)

  p_ym_0_noncentered_IFs <- matrix(NA, nrow = n, ncol = NROW(my_values) )
  p_ym_0_centered_IFs <- matrix(NA, nrow = n, ncol = NROW(my_values) )

  p_ym_1_noncentered_IFs <- matrix(NA, nrow = n, ncol = NROW(my_values) )
  p_ym_1_centered_IFs <- matrix(NA, nrow = n, ncol = NROW(my_values) )

  k <- length(mvalues)

  p_m_0_noncentered_IFs <- matrix(NA, nrow = n, ncol = k )
  p_m_0_centered_IFs <- matrix(NA, nrow = n, ncol = k )

  p_m_1_noncentered_IFs <- matrix(NA, nrow = n, ncol = k )
  p_m_1_centered_IFs <- matrix(NA, nrow = n, ncol = k )

# Randomised design (reg_formula is NULL)

  if (is.null(reg_formula)){

    for(i in 1:NROW(my_values)){
      p_ym_0_indicators <- (yvec == my_values$y[i]) &
        (mvec == my_values$m[i]) &
        (dvec == 0)

      p_ym_0_noncentered_IFs[,i] <- p_ym_0_indicators / (n0/n)
      p_ym_0_centered_IFs[,i] <- (dvec==0) * (p_ym_0_indicators - mean(p_ym_0_noncentered_IFs[,i])) / (n0/n)

      p_ym_1_indicators <- (yvec == my_values$y[i]) &
        (mvec == my_values$m[i]) &
        (dvec == 1)

      p_ym_1_noncentered_IFs[,i] <- p_ym_1_indicators / (n1/n)
      p_ym_1_centered_IFs[,i] <- (dvec==1) * (p_ym_1_indicators - mean(p_ym_1_noncentered_IFs[,i])) / (n1/n)
    }

    for(i in 1:k){
      p_m_0_indicators <-
        (mvec == mvalues[i]) &
        (dvec == 0)

      p_m_0_noncentered_IFs[,i] <- p_m_0_indicators / (n0/n)
      p_m_0_centered_IFs[,i] <- (dvec==0) * (p_m_0_indicators - mean(p_m_0_noncentered_IFs[,i])) / (n0/n)

      p_m_1_indicators <-
        (mvec == mvalues[i]) &
        (dvec == 1)

      p_m_1_noncentered_IFs[,i] <- p_m_1_indicators / (n1/n)
      p_m_1_centered_IFs[,i] <- (dvec==1) * (p_m_1_indicators - mean(p_m_1_noncentered_IFs[,i])) / (n1/n)
    }


  } else{
    # Non_experimental design
    df0 <- df
    df0[[d]] <- 0    # all units set to control
    df1 <- df
    df1[[d]] <- 1    # all units set to treated

    for(idx in 1:NROW(my_values)){
      df$lhs <- as.numeric(yvec == my_values$y[idx] & mvec == my_values$m[idx])

      df$lhs0 <- (df[[d]] - 1) * df$lhs
      df$lhs1 <- df[[d]] * df$lhs

      if(var(df$lhs0) == 0){
        p_ym_0_centered_IFs[, idx] <- 0
      }else{
        fml0 <- stats::as.formula( base::paste("lhs0", reg_formula) )
        reg0 <- fixest::feols(fml = fml0, data = df)

        S <- sandwich::estfun(reg0)  # n × p
        B <- sandwich::bread(reg0)   # the “bread”; for OLS it corresponds to (X'X)^(-1) up to n-scaling
        n <- nrow(S)

        # Observation-level influence function (n × p)
        IF <- (S %*% t(B))

        if(d %in% colnames(IF)){
          p_ym_0_centered_IFs[, idx] <- IF[, d]
        }else if( base::max( base::grepl(pattern = base::paste0(d," ="), x = colnames(IF)) ) >0 ){
          p_ym_0_centered_IFs[, idx] <- IF[, base::which(base::grepl(pattern = base::paste0(d," ="), x = colnames(IF)))]
        }

      }


        if(var(df$lhs1) == 0){
          p_ym_1_centered_IFs[, idx] <- 0
        }else{
          fml1 <- stats::as.formula( base::paste("lhs1", reg_formula) )
          reg1 <- fixest::feols(fml = fml1, data = df)

          S <- sandwich::estfun(reg1)  # n × p
          B <- sandwich::bread(reg1)   # the “bread”; for OLS it corresponds to (X'X)^(-1) up to n-scaling
          n <- nrow(S)

          # Observation-level influence function (n × p)
          IF <- (S %*% t(B))

          if(d %in% colnames(IF)){
            p_ym_1_centered_IFs[, idx] <- IF[, d]
          }else if( base::max( base::grepl(pattern = base::paste0(d," ="), x = colnames(IF)) ) >0 ){
            p_ym_1_centered_IFs[, idx] <- IF[, base::which(base::grepl(pattern = base::paste0(d," ="), x = colnames(IF)))]
          }

      }
    }


    for (i in 1:k){
      idx <- which(my_values$m == mvalues[i])
      if(length(idx)){
        p_m_0_centered_IFs[,i] <- rowSums(p_ym_0_centered_IFs[,idx,drop=FALSE])
        p_m_1_centered_IFs[,i] <- rowSums(p_ym_1_centered_IFs[,idx,drop=FALSE])
        p_m_0_noncentered_IFs[,i] <- rowSums(p_ym_0_noncentered_IFs[,idx,drop=FALSE])
        p_m_1_noncentered_IFs[,i] <- rowSums(p_ym_1_noncentered_IFs[,idx,drop=FALSE])
      }
    }
  }


  if (inequalities_only) {


    beta.obs_centered_IFs <- cbind(
      p_m_0_centered_IFs,
      p_m_1_centered_IFs,
      -p_m_0_centered_IFs,
      -p_m_1_centered_IFs,
      p_ym_1_centered_IFs - p_ym_0_centered_IFs)
  } else {

    beta.obs_centered_IFs <- cbind(
      p_m_0_centered_IFs,
      p_m_1_centered_IFs,
      p_ym_1_centered_IFs - p_ym_0_centered_IFs)
  }

  if (exploit_binary_m) {
    num_yvals <- length(unique(yvec))
    return(list(beta.obs_centered_IFs =
                  cbind((p_ym_0_centered_IFs - p_ym_1_centered_IFs)[,1:num_yvals], (p_ym_1_centered_IFs - p_ym_0_centered_IFs)[,(num_yvals+1):(2*num_yvals)])))
  }

  return(list(
              beta.obs_centered_IFs = beta.obs_centered_IFs,
              p_ym_0_centered_IFs = p_ym_0_centered_IFs,
              p_ym_1_centered_IFs = p_ym_1_centered_IFs,
              p_m_0_centered_IFs = p_m_0_centered_IFs,
              p_m_1_centered_IFs = p_m_1_centered_IFs
  ))

}

analytic_variance <-
  function(yvec, dvec, mvec, df, d, reg_formula = NULL,
           clustervec = seq(from = 1, to = length(yvec)),
           my_values, mvalues = unique(my_values$m),
           inequalities_only = TRUE, exploit_binary_m = FALSE){

    IF_out <- get_IFs(yvec = yvec,
                   dvec = dvec,
                   mvec = mvec,
                   df = df,
                   d = d,
                   reg_formula = reg_formula,
                   my_values = my_values,
                   mvalues = mvalues,
                   inequalities_only = inequalities_only,
                   exploit_binary_m = exploit_binary_m)
    IFs <- IF_out$beta.obs_centered_IFs


    #Sum the IFs within cluster
    IFs_clustered <- base::rowsum(x = IFs,
                                  group = clustervec)

    n <- length(yvec)
    c <- length(unique(clustervec))

    vcv <- cov(IFs_clustered * c/n) / NROW(IFs_clustered)

    return(vcv)
  }
