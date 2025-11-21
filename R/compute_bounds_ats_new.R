

#' @description
#' Given estimates of the CDF of a random variable Y, this function computes the expectation
#' of Y conditional on it being truncated to the top (or bottom) frac fraction of its distribution
#' For example, if frac = 0.5, we compute F^-1(U) | U in [0.5,1], where F is the CDF of Y
#' @param ecdf_table A data.frame containing columns y and cdf, where cdf denotes the estimated CDF of Y at y
#' @param frac The fraction in the tail over which we're taking the distribution
#' @param upper If TRUE, we look at the upper tail; if FALSE, we look at the lower tail. If NULL, we return a list with lb and ub.
#' @param num_gridpoints The number of gridpoints to use for numerical calculation

trimmed_expectation_from_cdf <- function(ecdf_table,
                                         frac,
                                         upper = NULL,
                                         num_gridpoints = 1e5) {
  # ---- Preprocess CDF ----
  stopifnot(is.data.frame(ecdf_table),
            all(c("y","cdf") %in% names(ecdf_table)))

  et <- ecdf_table[order(ecdf_table$y), , drop = FALSE]
  et$cdf <- as.numeric(et$cdf)
  et$cdf[!is.finite(et$cdf)] <- 0
  et$cdf <- pmin(1, pmax(0, et$cdf))  # clamp to [0,1]
  if (nrow(et)) {
    et$cdf <- cummax(et$cdf)          # enforce nondecreasing
    if (max(et$cdf, na.rm = TRUE) > 0) {
      et$cdf[nrow(et)] <- 1           # if any mass, force last point to 1
    }
  }

  # robust inverse CDF: returns finite y for all u in [0,1]
  inv_cdf <- function(u) {
    if (!nrow(et)) return(NA_real_)
    u <- max(0, min(1, as.numeric(u)))  # clamp u

    cdf <- et$cdf
    y   <- et$y

    # if CDF never reaches u -> return largest y (avoid empty set)
    if (max(cdf, na.rm = TRUE) < u) return(max(y, na.rm = TRUE))
    # if already >= u at the start -> return smallest y
    if (cdf[1] >= u) return(y[1])

    idx <- which(cdf >= u)[1]
    if (is.na(idx)) max(y, na.rm = TRUE) else y[idx]
  }

  # compute both tails if requested
  if (is.null(upper)) {
    return(list(
      lb = trimmed_expectation_from_cdf(et, frac, upper = FALSE, num_gridpoints = num_gridpoints),
      ub = trimmed_expectation_from_cdf(et, frac, upper = TRUE,  num_gridpoints = num_gridpoints)
    ))
  }

  # ---- u-grid construction with numerical guards ----
  frac <- as.numeric(frac)
  if (!is.finite(frac)) frac <- 0
  frac <- max(0, min(1, frac))

  eps <- sqrt(.Machine$double.eps)      # tiny gap to avoid exactly 1
  if (upper) {
    u_start <- max(0, 1 - frac)
    u_end   <- 1 - eps
    if (u_start > u_end) u_start <- u_end
  } else {
    u_start <- 0
    u_end   <- min(frac, 1 - eps)
    if (u_end < u_start) u_end <- u_start
  }

  if (num_gridpoints <= 1) {
    u_grid <- u_start
  } else {
    u_grid <- seq(from = u_start, to = u_end, length.out = num_gridpoints)
  }

  mean(vapply(u_grid, inv_cdf, numeric(1)))
}



#' @title Bounds on direct effect for always-takers
#' @description Computes bounds on E[Y(1,k) - Y(1,k) | G = kk], the average direct effect for k-ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param at_group The AT group of interest. Default is 1, so we compute means for units with M(1)=M(0)=1.
#' @param max_defier_share Bound on the proportion of defiers in the population. Default is 0 which indicates that the monotonicity constraint is imposed.
#' @param allow_min_defiers If the bound on defiers (max_defier_share) is inconsistent with the data, proceed by allowing the
#'   minimum number of defiers compatible with the data. Otherwise, throw an error. Default is TRUE.
#' @param num_gridpoints (Optional.) The number of gridpoints used in evaluating the integral. Higher is more accurate but more computationally costly
#' @importFrom "stats" "quantile"
#' @export

compute_bounds_ats_new <- function(df,
                               d,
                               m,
                               y,
                               at_group = 1,
                               max_defier_share = 0,
                               allow_min_defiers = TRUE,
                               reg_formula = NULL,
                               num_gridpoints = 10^5){


  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]
  n <- nrow(df)


  # compute theta_{kk}^{min} using the SAME LP structure as in lb_frac_affected,
  # but with a different objective (minimize \theta_{kk} instead of \nu_k).
  # replace Y|M=1,D with Y|M=k,D
  # replace at_share_m1d1 with \theta_{kk,min} / P(M=k|D=1) and at_share_m1d0 with \theta_{kk}^{min} / P(M=k | D=0)

  yvalues <- sort(unique(yvec))

  # row equality that works for scalar or multi-col m
  .row_equals <- function(x,y){
    xv <- as.vector(as.matrix(x))
    yv <- as.vector(if (is.matrix(y)||is.data.frame(y)) as.matrix(y) else rep(y, length(xv)))
    all(xv == yv)
  }

  mdf <- df[,m,drop = FALSE]
  mvalues <- unique(mdf)
  if (!is.data.frame(mvalues)) mvalues <- as.data.frame(mvalues)

  #Create a data frame containing all combinations of the row numbers of mvalues (same as lb_frac_affected)
  mvalues_df <- base::expand.grid(1:NROW(mvalues), 1:NROW(mvalues))
  m0_types <- mvalues[mvalues_df[,1],, drop = FALSE]
  m1_types <- mvalues[mvalues_df[,2],, drop = FALSE]

  #Determine which rows of m0_types are weakly less than the corresponding row in m1_types
  monotonic_types <- base::sapply(1:NROW(m0_types),
                                  FUN = function(i){
                                    base::all(m0_types[i,] <= m1_types[i,])})

  defier_types <- 1-monotonic_types

  if (is.null(reg_formula)) {
    p_m_1_fn <- function(mvalue){
      mdf1 <- mdf[dvec == 1,, drop = FALSE]
      M_equals_mvalue <- base::sapply(1:NROW(mdf1), function(i) .row_equals(mdf1[i,], mvalue))
      mean(M_equals_mvalue)
    }
    p_m_0_fn <- function(mvalue){
      mdf0 <- mdf[dvec == 0, , drop = FALSE]
      M_equals_mvalue <- base::sapply(1:NROW(mdf0), function(i) .row_equals(mdf0[i,], mvalue))
      mean(M_equals_mvalue)
    }

    p_m_1 <- base::apply(mvalues, 1, p_m_1_fn)
    p_m_0 <- base::apply(mvalues, 1, p_m_0_fn)
  } else {
    p_m_1 <- p_m_0 <- rep(NA_real_, nrow(mvalues))
  }


  regression_probs <- NULL

  if (is.null(reg_formula)) {
    max_p_diffs_list <- compute_max_p_difference(
      dvec = dvec,
      mdf  = mdf,
      yvec = yvec,
      wvec = rep(1, length(yvec)),
      continuous_Y = FALSE
    )
  } else {

    max_p_diffs_list <- compute_max_p_difference_reg(
      dvec = dvec,
      mdf  = mdf,
      yvec = yvec,
      wvec = rep(1, length(yvec)),
      continuous_Y = FALSE,
      df   = df,
      d    = d,
      y    = y,
      reg_formula = reg_formula
    )
    regression_probs <- max_p_diffs_list
  }
  max_p_diffs <- max_p_diffs_list$max_p_diffs

  if (!is.null(regression_probs)) {
    my_values_reg <- regression_probs$my_values
    p_ym_0_vec <- regression_probs$p_ym_0_vec
    p_ym_1_vec <- regression_probs$p_ym_1_vec

    p_m_0 <- vapply(
      seq_len(nrow(mvalues)),
      function(idx) {
        sum(p_ym_0_vec[my_values_reg$m == idx])
      },
      numeric(1)
    )

    p_m_1 <- vapply(
      seq_len(nrow(mvalues)),
      function(idx) {
        sum(p_ym_1_vec[my_values_reg$m == idx])
      },
      numeric(1)
    )
  }


  # build the same constraint matrices as in lb_frac_affected
  m1_marginals_constraints_matrix <-
    Reduce(rbind,
           base::lapply(X = 1:NROW(mvalues),
                        FUN = function(i){base::apply(m1_types, 1, function(x){ base::all(x == mvalues[i,]) })}))

  #Add a matrix of zeros to m1_marignals_constraints_matrix with length NROW(mvalues)
  m1_marginals_constraints_matrix <- base::cbind(m1_marginals_constraints_matrix,
                                                 matrix(0,
                                                        nrow = NROW(m1_marginals_constraints_matrix),
                                                        ncol = NROW(mvalues)))

  ##Create a matrix corresponding to the constraint that the sum of the m1_types must equal the marginal distribution of M among D=0
  #Create a matrix where each row i is a NROW(m0_types) length vector where the jth element indicates if the jth row of m0_types equals the ith row of mvalues
  m0_marginals_constraints_matrix <-
    Reduce(rbind,
           base::lapply(X = 1:NROW(mvalues),
                        FUN = function(i){base::apply(m0_types, 1, function(x){ base::all(x == mvalues[i,]) })}))

  #Add a matrix of zeros to m0_marignals_constraints_matrix
  m0_marginals_constraints_matrix <- base::cbind(m0_marginals_constraints_matrix,
                                                 matrix(0,
                                                        nrow = NROW(m0_marginals_constraints_matrix),
                                                        ncol = NROW(mvalues)))

  #Create a matrix (really, vector) that bounds the total defiers share
  defiers_constraints_matrix <- c(defier_types, rep(0,NROW(mvalues)))

  ## We now check feasibility of the program by computing the minimal defier share
  # consistent with the constraints
  feasibility_lp <-
    Rglpk::Rglpk_solve_LP(obj = defiers_constraints_matrix, #obj is sum of defiers shares
                          mat = rbind(m1_marginals_constraints_matrix,
                                      m0_marginals_constraints_matrix),
                          rhs = c(p_m_1,p_m_0),
                          dir = rep("==", NROW(m1_marginals_constraints_matrix)*2),
                          max = FALSE
    )

  if(feasibility_lp$status == 1){
    warning("Error in checking feasibility. Proceed with caution")
  }else{

    min_defier_share <- feasibility_lp$optimum

    if(min_defier_share > max_defier_share){
      if(allow_min_defiers){
        max_defier_share <- min_defier_share + 10^(-6) #update max defier share to the optimum plus a small tolerance
        warning(paste0("The data is incompatible with the specified max_defiers_share.\n",
                       "                           Setting this to the min value compatible with the data:",
                       max_defier_share))
      }else{
        stop(paste0("The data is incompatible with the specified max_defiers_share. \n",
                    "                        The specified value is ",
                    max_defier_share, " but the min compatible with the data is ", min_defier_share))
      }
    }
  }

  ##Create a constraint corresponding to the constraint that max_p_diffs must be greater than or equal to the sum of all complier types that end up at a given mvalue
  #Create a matrix where each row i is a NROW(m1_types) length vector where the jth element indicates if the jth row of m1_types equals the ith row of mvalues and the jth row of m1_types is not exactly equal to the jth row of m0_types
  complier_constraints_matrix <-
    Reduce(rbind, base::lapply(1:NROW(mvalues), function(m_index) {
      base::sapply(1:NROW(m1_types), function(s) {
        row_equals(m1_types[s, ], mvalues[m_index, ]) &&
          !row_equals(m0_types[s, ], mvalues[m_index, ])
      })
    }))

  complier_constraints_matrix <- base::cbind(
    complier_constraints_matrix,
    diag(NROW(mvalues))
  )

  #Combine the constraint matrices
  constraints_matrix <- base::rbind(m1_marginals_constraints_matrix,
                                    m0_marginals_constraints_matrix,
                                    complier_constraints_matrix,
                                    defiers_constraints_matrix)

  #Combine the constants associated with the matrices
  rhs_vec <- c(p_m_1, p_m_0, max_p_diffs, max_defier_share)

  #Specify the direction of the equalities/inequalities
  dir <- c(rep("==", 2*NROW(m1_marginals_constraints_matrix)),
           rep(">=", NROW(max_p_diffs)),
           "<=")

  ##Create a function that returns the objective vectors for a single at_group
  ##Our objective is theta_kk
  compute_obj_vectors <- function(at_group){


    #Find the index i such that the ith row of m0_types and m1_types are equal to at_group
    at_group_index <- which(purrr::map_lgl(.x = 1:NROW(m0_types),
                                           .f = ~row_equals(m0_types[.x,], at_group) &
                                             row_equals(m1_types[.x,], at_group)))

    #Thus denominator is a basis vector with one in pos corresponding to thetakk
    obj_denominator <- rep(0, NCOL(constraints_matrix))
    obj_denominator[at_group_index] <- 1

    #Find the index i such that the ith row of mvalues equals at_group
    at_group_index <- which(purrr::map_lgl(.x = 1:NROW(mvalues),
                                           .f = ~row_equals(mvalues[.x,],at_group)))

    #Numerator is a basis vector with one in pos corresponding to thetakk TVk
    obj_numerator <- rep(0, NCOL(constraints_matrix))
    obj_numerator[NROW(m1_types) + at_group_index] <- 1

    return(list(obj_numerator = obj_numerator,
                obj_denominator = obj_denominator))
  }

  #for at_group = k
  obj_denominator <- compute_obj_vectors(at_group = at_group)$obj_denominator
  # obj_numerator is not needed here, we don't solve the fractional LP in this function

  # LP to get \theta_{kk}^{min} with same constraints
  theta_lp <- Rglpk::Rglpk_solve_LP(
    obj    = obj_denominator,
    mat    = constraints_matrix,
    dir    = dir,
    rhs    = rhs_vec,
    bounds = NULL,
    max    = FALSE
  )

  theta_kk_min <- as.numeric(theta_lp$optimum)


  at_group_index <- which(purrr::map_lgl(
    .x = 1:NROW(mvalues),
    .f = ~.row_equals(mvalues[.x, , drop = FALSE], at_group)
  ))

  if (!length(at_group_index)) {
    warning("The requested at_group does not appear in the mediator support; returning NA bounds.")
    return(data.frame(lb = NA_real_, ub = NA_real_))
  }

  # Replace shares: θ_{kk}^{min} / P(M=k | D=d)
  if (is.null(reg_formula)) {
    p_mk_d1 <- mean(mvec[dvec == 1] == at_group)
    p_mk_d0 <- mean(mvec[dvec == 0] == at_group)
  } else {
    p_mk_d1 <- p_m_1[at_group_index]
    p_mk_d0 <- p_m_0[at_group_index]
  }
  if (p_mk_d1 <= 0 || p_mk_d0 <= 0) {
    warning("P(M=at_group | D=d) is zero for some d; returning NA bounds.")
    return(data.frame(lb = NA_real_, ub = NA_real_))
  }


  # trimming shares per Lemma A.3
  at_share_mkd1 <- theta_kk_min / p_mk_d1
  at_share_mkd0 <- theta_kk_min / p_mk_d0

  # ---- Build CDFs for Y | M=at_group, D=d ----
  # If reg_formula is NULL -> use ecdfs
  # If reg_formula is set -> build regression-inferred CDFs (discrete Y only).
  if(is.null(reg_formula)){

  ecdf_y_mkd1_fn <- stats::ecdf(yvec[ mvec == at_group & dvec ==1 ]) #this returns a fn y -> CDF
  ecdf_y_mkd1 <- ecdf_y_mkd1_fn(yvalues) #this is a vector of CDFs evaluated at yvalues

  #Similarly we compute the CDF of Y | M=k,D=0
  ## In this case, this trimming argument is trivial, but setting it up this week for generalization
  ecdf_y_mkd0_fn <- stats::ecdf(yvec[ mvec == at_group & dvec ==0 ]) #this returns a fn y -> CDF
  ecdf_y_mkd0 <- ecdf_y_mkd0_fn(yvalues) #this is a vector of CDFs evaluated at yvalues
  } else{

    if (!length(at_group_index)) {
      warning("The requested at_group does not appear in the mediator support; returning NA bounds.")
      return(data.frame(lb = NA_real_, ub = NA_real_))
    }

    m_match <- apply(mdf, 1, function(row) .row_equals(row, at_group))

    extract_joint <- function(prob_vec) {
      vapply(
        yvalues,
        function(y_val) {
          idx <- which(my_values_reg$y == y_val & my_values_reg$m == at_group_index)
          if (length(idx) == 0) 0 else prob_vec[idx]
        },
        numeric(1)
      )
    }

    joint1 <- pmax(0, extract_joint(p_ym_1_vec))
    joint0 <- pmax(0, extract_joint(p_ym_0_vec))

    build_conditional_pmf <- function(joint, y_cell) {
      s <- sum(joint)
      if (!is.finite(s) || s <= .Machine$double.eps) {
        if (length(y_cell) == 0L) return(rep(0, length(yvalues)))
        tab <- table(factor(y_cell, levels = yvalues))
        return(as.numeric(tab) / sum(tab))
      }
      pmf <- joint / s
      pmf[pmf < 0] <- 0
      s2 <- sum(pmf)
      if (s2 > 0) pmf <- pmf / s2
      pmf
    }

    pmf1 <- build_conditional_pmf(joint1, y_cell = yvec[m_match & dvec == 1])
    pmf0 <- build_conditional_pmf(joint0, y_cell = yvec[m_match & dvec == 0])

    validate_cdf <- function(cdfv) {
      if (!length(cdfv)) return(cdfv)
      cdfv[!is.finite(cdfv)] <- 0
      cdfv <- pmin(1, pmax(0, cdfv))
      cdfv <- cummax(cdfv)
      cdfv[length(cdfv)] <- 1
      cdfv
    }

    ecdf_y_mkd1 <- validate_cdf(cumsum(pmf1))
    ecdf_y_mkd0 <- validate_cdf(cumsum(pmf0))
    }



  # The ub is the top at_share_mkd1 of the Y|M=k,D=1 distribution
  # The lb is the bottom at_share_mkd1 of the Y|M=k,D=1 distribution
  bounds_y1 <- trimmed_expectation_from_cdf(ecdf_table = data.frame(y=yvalues, cdf = ecdf_y_mkd1),
                                            frac = at_share_mkd1,
                                            upper = NULL,
                                            num_gridpoints = num_gridpoints)

  # The ub is the top at_share_mkd0 of the Y|M=k,D=1 distribution
  # The lb is the bottom at_share_mkd0 of the Y|M=k,D=1 distribution
  bounds_y0 <- trimmed_expectation_from_cdf(ecdf_table = data.frame(y=yvalues, cdf = ecdf_y_mkd0),
                                            frac = at_share_mkd0,
                                            upper = NULL,
                                            num_gridpoints = num_gridpoints)

  ##Now compute bounds on ADE for always-takers
  # LB is lb for y1 minus ub for y0
  lb <- bounds_y1$lb - bounds_y0$ub
  ub <- bounds_y1$ub - bounds_y0$lb

  return(data.frame(lb=lb, ub=ub))
}
