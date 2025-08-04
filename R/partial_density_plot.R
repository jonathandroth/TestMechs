#' @title Create plots of partial densities
#' @description Plots f_{Y,M=m | D=1} and f_{Y,M=m | D=0} for m equal to 0 or 1
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param numGridPoints Number of points used in grid for graphing
#' @param plot_nts (Optional) If TRUE, we plot f_{Y,M=0 | D=1} and f_{Y,M=0 | D=0}. Otherwise, we plot f_{Y,M=1 | D=1} and f_{Y,M=1 | D=0}. Default is FALSE
#' @param density_1_label (Optional) The label on the plot for the d=1 density.
#' @param density_0_label (Optional) The label on the plot for the d=0 density.
#' @param continuous_Y (Optional) Should Y be treated as continuous for density estimation. Default is TRUE. Use FALSE for discrete Y
#' @param reg_formula (Optional) Regression formula for observational adjustment  
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @return A ggplot object showing partial densities
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_bar

#' @export

partial_density_plot <- function(df,
                                 d,
                                 m,
                                 y,
                                 numGridPoints = 10000,
                                 plot_nts = FALSE,
                                 density_1_label = "f(Y,M=1|D=1)",
                                 density_0_label = "f(Y,M=1|D=0)",
                                 num_Ybins = NULL,
                                 reg_formula = NULL,
                                 continuous_Y = base::ifelse(is.null(num_Ybins),
                                                             TRUE,FALSE)){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  if (length(unique(df[[m]])) > 2) stop("This plot should be used only when M is binary")
  #If plot_nts = TRUE, re-run with m -> 1-m, d -> 1-d and flip the labels
  if(plot_nts == TRUE){
    df[[m]] <- 1- df[[m]]
    df[[d]] <- 1- df[[d]]

    # If the default labels are given, re-do the defaults
    # Note that D=1 after flipping corresponds to the original D=0
    # and likewise for M
    if(density_1_label == "f(Y,M=1|D=1)" & density_0_label == "f(Y,M=1|D=0)"){
      density_1_label <- "f(Y,M=0|D=1)"
      density_0_label <- "f(Y,M=0|D=0)"
    }else{
      #If custom labels are given, flip which one corresponds to D=1 and D=0
      density_1_label_old <- density_1_label
      density_1_label <- density_0_label
      density_0_label <- density_1_label_old
    }

    return(
      partial_density_plot(df = df,
                           d = d,
                           m = m,
                           y = y,
                           numGridPoints = numGridPoints,
                           plot_nts = FALSE,
                           density_1_label = density_1_label,
                           density_0_label = density_0_label,
                           num_Ybins = num_Ybins,
                           reg_formula = reg_formula,
                           continuous_Y = continuous_Y))
  }

  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
df[[y]] <- yvec
  }

  partial_densities_and_shares <- compute_partial_densities_and_shares(df = df,
                                                                       d = d,
                                                                       m = m,
                                                                       y = y,
                                                                       n = numGridPoints,
                                                                       continuous_Y = continuous_Y,
                                                                       reg_formula = reg_formula)


  if(continuous_Y == TRUE){
  f_partial11 <- partial_densities_and_shares$f_partial11
  f_partial01 <- partial_densities_and_shares$f_partial01

  ygrid <- seq(from = base::min(yvec) - 1* stats::sd(yvec),
               to = base::max(yvec) + 1* stats::sd(yvec) ,
               length.out = numGridPoints)

  partial11_grid <- f_partial11(ygrid)
  partial01_grid <- f_partial01(ygrid)

  partial_density_df <-
  dplyr::bind_rows(
    data.frame(y = ygrid, dens = partial11_grid, Partial.Density = density_1_label),
    data.frame(y = ygrid, dens = partial01_grid, Partial.Density = density_0_label)
    )

  partial_density_plot <-
  partial_density_df %>%
    ggplot(aes(x = y,
               y = dens,
               color = Partial.Density,
               linetype = Partial.Density)) +
    geom_line() +
    xlab("Y") +
    ylab("Partial Density")
  }else{
    pmf_partial11 <- partial_densities_and_shares$pmf_partial11
    pmf_partial01 <- partial_densities_and_shares$pmf_partial01
    yvalues <- partial_densities_and_shares$yvalues

    partial_density_df <-
      dplyr::bind_rows(
        data.frame(y = yvalues, dens = pmf_partial11, Partial.Density = density_1_label),
        data.frame(y = yvalues, dens = pmf_partial01, Partial.Density = density_0_label)
      )


    partial_density_plot <-
      partial_density_df %>%
      ggplot(aes(x = y,
                 y = dens,
                 color = Partial.Density,
                 fill= Partial.Density,
                 linetype = Partial.Density)) +
      geom_bar(stat="identity", position = "dodge", width = 0.7) +
      xlab("Y") +
      ylab("Partial PMF")


  }

  return(partial_density_plot)
}

compute_partial_densities_and_shares <- function(df, 
                                                 d, 
                                                 m, 
                                                 y, 
                                                 w = NULL,
                                                 continuous_Y=TRUE,
                                                 reg_formula = NULL,
                                                 ...){
    yvec <- df[[y]]
    dvec <- df[[d]]
    mvec <- df[[m]]
    
    if(is.null(w)){
      wvec <- rep(1, NROW(df))
    }else{
      wvec <- df[[w]]
    }
    
    # Randomized branch
    if(is.null(reg_formula)){
      
    frac_compliers <- stats::weighted.mean( mvec[dvec == 1], w = wvec[dvec == 1] ) -
      stats::weighted.mean( mvec[dvec == 0], w = wvec[dvec == 0])
    frac_ats <- stats::weighted.mean( mvec[dvec == 0], w= wvec[dvec == 0] )
    theta_ats <- frac_ats / (frac_compliers + frac_ats) #fraction among Cs/ATs
    
    ats_untreated_index <- (dvec == 0) & (mvec == 1) #these are ATs when untreated
    ats_treated_index <- (dvec == 1) & (mvec == 1) #these are ATs or Cs
    
    y_ats_treated <- yvec[ats_treated_index]
    y_ats_untreated <- yvec[ats_untreated_index]
    
    w_ats_treated <- wvec[ats_treated_index]
    w_ats_untreated <- wvec[ats_untreated_index]
    
    #The density function doesn't normalize weights, so normalize these
    w_ats_treated <- w_ats_treated/sum(w_ats_treated)
    w_ats_untreated <- w_ats_untreated/sum(w_ats_untreated)
    
    if(continuous_Y){
      dens_y_ats_treated <- get_density_fn(x = y_ats_treated, weights = w_ats_treated, ...)
      dens_y_ats_untreated <- get_density_fn(x = y_ats_untreated, weights = w_ats_untreated, ...)
      
      f_partial11 <- function(y){ (frac_ats + frac_compliers) * dens_y_ats_treated(y) }
      f_partial01 <- function(y){ frac_ats  * dens_y_ats_untreated(y) }
      
      resultsList <-
        list(frac_compliers = frac_compliers,
             frac_ats = frac_ats,
             theta_ats = theta_ats,
             f_partial11 = f_partial11,
             f_partial01 = f_partial01)
      
      return(resultsList)
    }else{
      yvalues <- unique(yvec)
      pmf_y_ats_treated <- purrr::map_dbl(.x = 1:length(yvalues),
                                          .f = ~stats::weighted.mean(x = y_ats_treated == yvalues[.x],
                                                                     w = w_ats_treated))
      
      pmf_y_ats_untreated <- purrr::map_dbl(.x = 1:length(yvalues),
                                            .f = ~stats::weighted.mean(x = y_ats_untreated == yvalues[.x],
                                                                       w = w_ats_untreated))
      
      pmf_partial_11 <- (frac_ats + frac_compliers) * pmf_y_ats_treated
      pmf_partial_01 <- (frac_ats) * pmf_y_ats_untreated
      
      resultsList <-
        list(frac_compliers = frac_compliers,
             frac_ats = frac_ats,
             theta_ats = theta_ats,
             pmf_partial11 = pmf_partial_11,
             pmf_partial01 = pmf_partial_01,
             yvalues = yvalues)
    }
    } else{
      # Regression branch
      if (continuous_Y){
        stop("reg_formula is only allowed when Y is discrete (set continuous_Y = FALSE or supply num_Ybins).")
      }else{
     
      # ---- ensure treatment on RHS ---------------------------------
      reg_formula_chr <- paste(deparse(reg_formula), collapse = " ")
      if (!grepl("\\btreat\\b", reg_formula_chr)) {
        warning("The treatment variable '", d,
                "' was not found in the provided reg_formula; ",
                "I have added it as a regressor. Please edit reg_formula if that was not your intention.")
        reg_formula_chr <- paste("~", d, "+", sub("^~", "", reg_formula_chr))
        reg_formula     <- as.formula(reg_formula_chr)
      }
      
      # ---- predicted P(M=1 | D=d)  ----------------------------------
      p_m1 <- .predict_mean(lhs_vec = as.integer(mvec == 1),
                            df      = df,
                            d_var   = d,
                            d_val   = 1,
                            reg_formula = reg_formula,
                            wvec    = wvec)
      p_m0 <- .predict_mean(lhs_vec = as.integer(mvec == 1),
                            df      = df,
                            d_var   = d,
                            d_val   = 0,
                            reg_formula = reg_formula,
                            wvec    = wvec)
      
      frac_compliers <- p_m1 - p_m0
      frac_ats       <- p_m0
      theta_ats      <- frac_ats / (frac_compliers + frac_ats)
      
      # ---- joint probabilities for each y ---------------------------
      yvalues <- sort(unique(yvec))
      
      pmf_partial11 <- purrr::map_dbl(
        yvalues,
        ~.predict_mean(lhs_vec = as.integer((yvec == .x) & (mvec == 1)),
                       df      = df,
                       d_var   = d,
                       d_val   = 1,
                       reg_formula = reg_formula,
                       wvec    = wvec)
      )
      
      pmf_partial01 <- purrr::map_dbl(
        yvalues,
        ~.predict_mean(lhs_vec = as.integer((yvec == .x) & (mvec == 1)),
                       df      = df,
                       d_var   = d,
                       d_val   = 0,
                       reg_formula = reg_formula,
                       wvec    = wvec)
      )
      
      # Return (discrete-Y only)
      resultsList <- list(frac_compliers = frac_compliers,
           frac_ats       = frac_ats,
           theta_ats      = theta_ats,
           pmf_partial11  = pmf_partial11,
           pmf_partial01  = pmf_partial01,
           yvalues        = yvalues)
    }
    }
    return(resultsList)
}

.predict_mean <- function(lhs_vec,
                          df,
                          d_var,
                          d_val,
                          reg_formula,
                          cluster = NULL,
                          wvec = NULL){
  
  iv_spec <- extract_iv(reg_formula, d_var)
  
  # add lhs to df
  df$lhs <- lhs_vec
  
  # variables needed for NA-drop
  need_vars <- unique(c(d_var, "lhs", iv_spec$controls, iv_spec$instr))
  need_vars <- need_vars[nzchar(need_vars)]
  if (!is.null(cluster)) need_vars <- unique(c(need_vars, cluster))
  
  keep_idx <- stats::complete.cases(df[, need_vars, drop = FALSE])
  if (!any(keep_idx)) return(NA_real_)
  
  df_clean <- df[keep_idx, , drop = FALSE]
  
  # if lhs is constant after NA-drop, just return that constant
  lhs_clean <- df_clean$lhs
  lhs_clean <- lhs_clean[!is.na(lhs_clean)]
  if (length(lhs_clean) == 0L || length(unique(lhs_clean)) == 1L) {
    const_mean <- mean(lhs_clean)  # 0 or 1 (or NA if empty)
    return(const_mean)
  }
  
  # build fixest formula
  if (iv_spec$is_iv) {
    ctrl  <- paste(iv_spec$controls, collapse = " + ")
    ctrl_rhs <- if (nzchar(ctrl)) ctrl else "1"
    endo <- paste(iv_spec$treat, collapse = " + ")
    instr <- paste(iv_spec$instr, collapse = " + ")
    fml <- as.formula(sprintf("lhs ~ %s | %s ~ %s", ctrl_rhs, endo, instr))
  } else {
    rhs <- paste(c(iv_spec$treat, iv_spec$controls), collapse = " + ")
    rhs <- if (nzchar(rhs)) rhs else "1"
    fml <- as.formula(paste("lhs ~", rhs))
  }
  
  # counterfactual D value
  df_cf <- df_clean
  df_cf[[d_var]] <- d_val
  
  # fit and predict (cluster only affects SEs, not needed here)
  reg <- fixest::feols(fml, data = df_clean)
  preds <- as.numeric(predict(reg, newdata = df_cf))
  
  # average predictions (optionally weighted) over the kept rows
  if (is.null(wvec)) {
    return(mean(preds, na.rm = TRUE))
  } else {
    w_clean <- wvec[keep_idx]
    # normalize to avoid scale affecting the mean
    sw <- sum(w_clean, na.rm = TRUE)
    if (isTRUE(sw == 0) || is.na(sw)) {
      return(mean(preds, na.rm = TRUE))
    } else {
      w_clean <- w_clean / sw
      return(stats::weighted.mean(preds, w_clean, na.rm = TRUE))
    }
  }
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

