

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
#' @param max_defiers_share Bound on the proportion of defiers in the population. Default is 0 which indicates that the monotonicity constraint is imposed.
#' @param num_gridpoints (Optional.) The number of gridpoints used in evaluating the integral. Higher is more accurate but more computationally costly
#' @importFrom "stats" "quantile"
#' @export

compute_bounds_ats_new <- function(df,
                               d,
                               m,
                               y,
                               at_group = 1,
                               max_defier_share = 0,
                               reg_formula = NULL,
                               num_gridpoints = 10^5){

  #features not yet implemented. Chen, you can fill these in and remove these errors
  # if(at_group != 1){stop("function is currently only implemented for at_group = 1")}
  # if(!is.null(reg_formula)){stop("reg_formula not yet implemented")}
  # if(!is.null(max_defier_share)){stop("function is currently only implemented when max_defier_share = NULL") }
  

  # If at_group = 0 (Never-takers) 
  # Flip D and M to convert never-takers into always-takers in a flipped scenario
  if (at_group == 0) {
    df[[m]] <- 1 - df[[m]]
    df[[d]] <- 1 - df[[d]]
    
    # Compute bounds in the flipped scenario treating flipped group as always-takers
    bounds_flipped <- compute_bounds_ats_new(df = df,
                                             d = d, 
                                             m = m, 
                                             y = y, 
                                             at_group = 1, 
                                             max_defier_share = max_defier_share, 
                                             reg_formula = reg_formula, 
                                             num_gridpoints = num_gridpoints)
    # Flip the sign of the bounds (because we flipped D)
    return(data.frame(
      lb = -bounds_flipped$ub,
      ub = -bounds_flipped$lb
    ))
  }

  
  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]
  n <- nrow(df)


  # feasibility check when m is a 0/1 indicator
  is_binary_vec <- function(v){
    u <- sort(unique(v[is.finite(v)]))
    length(u) <= 2 && all(u %in% c(0,1))
  }
  
  # Try to infer if m == 1{ parent == some value } for an ordinal/numeric parent column.
  infer_parent_of_indicator <- function(df, mname){
    mv <- df[[mname]]
    if (!is_binary_vec(mv)) return(NULL)
    
    # We look for a single column 'parent' and value 'v' such that (parent==v) <-> (m==1)
    for (cn in names(df)) {
      if (cn == mname) next
      col <- df[[cn]]
      # only consider numeric/integer/ordered-factor columns as plausible parents
      if (!(is.numeric(col) || is.integer(col) || is.ordered(col))) next
      
      pos <- which(mv == 1 & is.finite(mv))
      if (!length(pos)) next
      uvals <- unique(col[pos])
      uvals <- uvals[is.finite(uvals)]
      if (length(uvals) != 1) next
      v <- uvals[1]
      
      # Check equivalence on rows where both sides are observed
      both_ok <- is.finite(mv) & is.finite(col)
      equiv <- all((col[both_ok] == v) == (mv[both_ok] == 1))
      if (!equiv) next
      
      # Determine if v is min or max value in the parent column
      parent_vals <- sort(unique(col[is.finite(col)]))
      is_min <- (v == min(parent_vals))
      is_max <- (v == max(parent_vals))
      return(list(parent = cn, v = v, is_min = is_min, is_max = is_max))
    }
    NULL
  }
  
  # If m is binary, try to detect if it's a middle-category indicator of some parent column.
  if (is_binary_vec(mvec)) {
    inf <- infer_parent_of_indicator(df, m)
    if (!is.null(inf)) {
      if (!(inf$is_min || inf$is_max)) {
        warning(sprintf(
          "Mediator '%s' looks like 1{%s == %s}, which is a middle category of '%s'. ",
          m, inf$parent, as.character(inf$v), inf$parent
        ),
        "Binary monotonicity is not valid for a middle-category dummy. ",
        "Please use the multi-valued mediator ('", inf$parent, "') with at_group = ", as.character(inf$v), ".")
        return(data.frame(lb = NA_real_, ub = NA_real_))
      }
      # else: it's min or max -> ok to proceed down the binary path
    }
    # If we couldn't infer a parent, we proceed (cannot reliably tell; keep previous behavior)
  }

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
  
  #Compute the marginal distribution of M among D=1
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
  

  
  # max_p_diffs
  
  # decide if Y is discrete:
  # - if reg_formula is provided, we REQUIRE discrete Y for the regression CDF path
  is_discrete_Y <- !( n / length(unique(df[[y]])) <= 30 )
  
  
  if (is.null(reg_formula)) {
    max_p_diffs_list <- compute_max_p_difference(
      dvec = dvec,
      mdf  = mdf,
      yvec = yvec,
      wvec = rep(1, length(yvec)),
      continuous_Y = !is_discrete_Y
    )
  } else {
    # regression-adjusted version (discrete Y only)
    if (!is_discrete_Y) {
      stop("reg_formula is supported only for discrete Y in compute_bounds_ats_new().")
    }
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
  }
  max_p_diffs <- max_p_diffs_list$max_p_diffs
  
  
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
  rhs_vec <- c(p_m_1, p_m_0, max_p_diffs, max_defiers_share)
  
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
  
  
  # Replace shares: θ_{kk}^{min} / P(M=k | D=d)
  p_mk_d1 <- mean(mvec[dvec == 1] == at_group)
  p_mk_d0 <- mean(mvec[dvec == 0] == at_group)
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
    # Regression-inferred CDFs (discrete Y only)
    yuniq <- sort(unique(yvec))
    if (!is_discrete_Y) {
      stop("reg_formula currently supported for discrete Y; Y variable might be continuous.")
    }
    
    # Match mediator pattern: M == at_group
    m_match <- apply(mdf, 1, function(row) row_equals(row, at_group))
    
    # Parse OLS/IV spec
    iv_spec <- extract_iv(reg_formula, d)
    
    # Allocate joint probabilities for each y and each d
    partial_pmf_d1 <- numeric(length(yvalues))
    partial_pmf_d0 <- numeric(length(yvalues))
    
    for (j in seq_along(yvalues)) {
      y_val <- yvalues[j]
      df$lhs <- as.integer(yvec == y_val & m_match)  # indicator {Y=y_val, M=at_group}
      
      # Build fixest formula
      if (iv_spec$is_iv) {
        ctrl  <- paste(iv_spec$controls, collapse = " + ")
        instr <- paste(iv_spec$instr,     collapse = " + ")
        fml   <- as.formula(sprintf("lhs ~ %s | %s ~ %s",
                                    ctrl,
                                    paste(iv_spec$treat, collapse = "+"),
                                    instr))
      } else {
        rhs <- paste(c(iv_spec$treat, iv_spec$controls), collapse = " + ")
        fml <- as.formula(paste("lhs ~", rhs))
      }
      
      # Vars, drop NA
      req <- c("lhs", d, iv_spec$controls, if (iv_spec$is_iv) iv_spec$instr else NULL)
      req <- unique(req[nzchar(req)])
      df_fit <- tidyr::drop_na(df[, req, drop = FALSE])
      if (nrow(df_fit) == 0) { partial_pmf_d1[j] <- 0; partial_pmf_d0[j] <- 0; next }
      
      # Fit model
      mod <- fixest::feols(fml, data = df_fit)
      
      # Counterfactual sets for D=1 and D=0
      new1 <- df_fit
      new1[[d]] <- 1
      new0 <- df_fit
      new0[[d]] <- 0
      
      pred1 <- as.numeric(predict(mod, newdata = new1))
      pred0 <- as.numeric(predict(mod, newdata = new0))
      
      # Average predicted joint probs: P(Y=y, M=k | D=d)
      partial_pmf_d1[j] <- mean(pred1, na.rm = TRUE)
      partial_pmf_d0[j] <- mean(pred0, na.rm = TRUE)
    }
    
    # Convert to conditional PMFs of Y | M=k, D=d
    # Nonnegative joint predictions (from your regression loop):
    joint1 <- pmax(0, as.numeric(partial_pmf_d1))  # ≈ P(Y=y, M=at_group | D=1)
    joint0 <- pmax(0, as.numeric(partial_pmf_d0))  # ≈ P(Y=y, M=at_group | D=0)
    
    # Build a VALID conditional PMF with safe empirical fallback
    build_conditional_pmf <- function(joint, y_cell) {
      s <- sum(joint)
      if (!is.finite(s) || s <= .Machine$double.eps) {
        # fallback: empirical conditional PMF of Y | M=at_group, D=d
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
    
    pmf1 <- build_conditional_pmf(joint1, y_cell = yvec[mvec == at_group & dvec == 1])
    pmf0 <- build_conditional_pmf(joint0, y_cell = yvec[mvec == at_group & dvec == 0])
    
    # Raw CDF vectors
    ecdf_y_mkd1 <- cumsum(pmf1)
    ecdf_y_mkd0 <- cumsum(pmf0)
    
    # --- Validator that takes ONLY the CDF vector and returns a fixed CDF vector ---
    .validate_cdf <- function(cdfv) {
      if (!length(cdfv)) return(cdfv)
      cdfv[!is.finite(cdfv)] <- 0
      cdfv <- pmin(1, pmax(0, cdfv))  # clamp to [0,1]
      cdfv <- cummax(cdfv)            # enforce nondecreasing
      cdfv[length(cdfv)] <- 1         # FORCE last point to 1
      cdfv
    }
    
    # Validate the CDF vectors
    ecdf_y_mkd1 <- .validate_cdf(ecdf_y_mkd1)
    ecdf_y_mkd0 <- .validate_cdf(ecdf_y_mkd0)
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