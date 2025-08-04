#' @title Bounds on fraction of always-takers affected
#'@description Computes bounds on the TV distance btwn Y(1,1) and Y(1,0) for always takers
#' This is a lower bound on the fraction of ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param continuous_Y (Optional) Whether Y should be treated as continuous, in which case kernel density is used, or discrete. Default is TRUE.
#@param method Either "density", to use a kernel density to compute TV, or "bins", to use a discrete approximation
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @export

compute_tv_ats <- function(df, d, m, y, w = NULL,
                           continuous_Y = base::ifelse(is.null(num_Ybins),TRUE,FALSE),
                           num_Ybins = NULL){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y,
                               w = w)


  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }

  dvec <- df[[d]]
  mvec <- df[[m]]

  if(is.null(w)){
    wvec <- rep(1, NROW(df))
  }else{
    wvec <- df[[w]]
  }

  #Use the general function to compute \int (partial_dens_1 - partial_dens0)_+
  max_p_diffs <-
  compute_max_p_difference(dvec = dvec,
                           mdf = data.frame(m= mvec),
                           yvec = yvec,
                           wvec = wvec,
                           continuous_Y = continuous_Y
                           )

  # Extract \int (partial_dens_1 - partial_dens_0)_+
  p_diff_1 <- max_p_diffs$max_p_diffs[which(max_p_diffs$mvalues == 1)]

  # What we want is 1/share_ats * (\int (partial_dens_1 - partial_dens_0)_+ - share_compliers)
  p_m_1 <- stats::weighted.mean(x = mvec[dvec == 1] == 1,
                                w = wvec[dvec == 1])
  p_m_0 <- stats::weighted.mean(x = mvec[dvec == 0] == 1,
                                w = wvec[dvec == 0])

  frac_compliers <- p_m_1 - p_m_0
  frac_ats <- p_m_0

  tv_lb <- 1/frac_ats * pmax(0,p_diff_1 - frac_compliers)

  # ### Below is our old, more direct computation that doesn't use functions from the multiple m case
  # partial_densities_and_shares <- compute_partial_densities_and_shares(df,d,m,y,w=w)
  #
  # f_partial01 <- partial_densities_and_shares$f_partial01
  # f_partial11 <- partial_densities_and_shares$f_partial11
  # theta_ats <- partial_densities_and_shares$theta_ats
  #
  # positive_part <- function(y){ base::pmax(y,0) }
  #
  # #Lower bound is 1/theta_{AT} \int positive_part( f_partial01 - f_partial11 )
  # #Implement this numerically over a grid of 10,000 points
  #
  # ygrid <- seq(from = base::min(yvec) - 3* stats::sd(yvec),
  #              to = base::max(yvec) + 3* stats::sd(yvec) ,
  #              length.out = 10000)
  #
  # TV_lb <-
  #   1/theta_ats *
  # base::sum( positive_part( f_partial01(ygrid) - f_partial11(ygrid) ) ) * base::diff(ygrid)[1]



  return(tv_lb)
}


#' @title Computes the breakdown number of defiers consistent with the sharp null
#' when M is binary
#'@description Computes the minimum fraction of the population that would need
#'to be defiers to be consistent with the sharp null.
#' Also computes the minimum ratio of defiers to compliers
#' @param df A data frame
#' @param d Name of the treatment variable in the df. Must be binary
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param continuous_Y (Optional) Whether Y should be treated as continuous, in which case kernel density is used, or discrete. Default is TRUE.
#@param method Either "density", to use a kernel density to compute TV, or "bins", to use a discrete approximation
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @export

compute_min_defiers_binaryM <- function(df, d, m, y, w = NULL,
                           continuous_Y = base::ifelse(is.null(num_Ybins),TRUE,FALSE),
                           num_Ybins = NULL){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y,
                               w = w)


  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }

  dvec <- df[[d]]
  mvec <- df[[m]]

  if(is.null(w)){
    wvec <- rep(1, NROW(df))
  }else{
    wvec <- df[[w]]
  }

  #Use the general function to compute \int (partial_dens_1 - partial_dens0)_+
  max_p_diffs <-
    compute_max_p_difference(dvec = dvec,
                             mdf = data.frame(m= mvec),
                             yvec = yvec,
                             wvec = wvec,
                             continuous_Y = continuous_Y
    )


  #Compute probabilities of M | D
  p_m_1 <- stats::weighted.mean(x = mvec[dvec == 1] == 1,
                                w = wvec[dvec == 1])
  p_m_0 <- stats::weighted.mean(x = mvec[dvec == 0] == 1,
                                w = wvec[dvec == 0])

  # Extract \int (partial_dens_{Y,M=1 | D=1} - partial_dens_{Y,M=1 | D=0})_+
  #Note that we actually want the negative part of this, but we have
  # (partial_dens_{Y,M=1 | D=1} - partial_dens_{Y,M=1 | D=0})_- =
  # (partial_dens_{Y,M=1 | D=1} - partial_dens_{Y,M=1 | D=0})_+ - (P(M=1 |D=1) - P(M=1|D=0))
  p_diff_1_reversed <- max_p_diffs$max_p_diffs[which(max_p_diffs$mvalues == 1)]
  p_diff_1 <- p_diff_1_reversed - (p_m_1 - p_m_0)

  # Extract \int (partial_dens_{Y,M=0 | D=1} - partial_dens_{Y,M=0 | D=0})_+
  p_diff_0 <- max_p_diffs$max_p_diffs[which(max_p_diffs$mvalues == 0)]

  #Minimum number of defiers is min of the p_diffs (or 0)
  d_lb <- pmax(0, p_diff_0, p_diff_1)

  #Number of compliers if P(M=1 | D=1) - P(M=1 | D=0) + (number of defiers)
  theta_c <- p_m_1 - p_m_0 + d_lb

  ratio <- d_lb / theta_c

  return(list(defier_lb = d_lb,
              ratio_lb = d_lb / theta_c))
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


get_density_fn <- function(x,...){

    dens <- stats::density(x=x,...)

    if("weights" %in% names(list(...))){
      warning("You are using weights with the continuous_Y option.
              R is weird about weighted densities, so this may not give you what you want.
              Try discretizing the outcome using num_Ybins instead.")
    }
    #Create a function that returns the density at any point
    # Return 0 if outside of range of dens$x
    # Otherwise, return the density at the largest pt in dens$x below y
    dens_fn <- function(y){
    index <- base::findInterval(y, dens$x)
    y_plus_zeros <- c(0,dens$y,0)
    return(y_plus_zeros[index+1])
  }

  return(dens_fn)
}

TV_distance_fn <- function(y1,
                           y2,
                           method = "density",
                           numbins = ceiling(length(y1)/20)){

  if(method == "bins"){
  binned_y <- dplyr::ntile( c(y1,y2), numbins)

  binned_y1 <- binned_y[(1:length(y1))]
  binned_y2 <- binned_y[-(1:length(y1))]

  y1_freq_table <- data.frame(y1 = binned_y1) %>%
                    dplyr::group_by(y1) %>%
                    dplyr::count() %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(dens1  = n/base::sum(n) )

    y2_freq_table <- data.frame(y2 = binned_y2) %>%
                    dplyr::group_by(y2) %>%
                    dplyr::count() %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(dens2  = n/sum(n) )

  joint_freq_table <- dplyr::full_join(y1_freq_table, y2_freq_table,
                                 by = c("y1" = "y2"))


  joint_freq_table <- joint_freq_table %>%
                        dplyr::mutate(dens1 = base::ifelse(is.na(dens1), 0 , dens1),
                               dens2 = base::ifelse(is.na(dens2), 0 , dens2),)

  TV <- 1/2 * base::sum(base::abs(joint_freq_table$dens1 - joint_freq_table$dens2))
  return(TV)
  }
  if(method == "density"){
    dens1 <- stats::density(x = y1)
    dens2 <- stats::density(x = y2)

    dens_fn <- function(dens,y){
      # if(y < min(dens$x) | y> max(dens$x)){
      #   #If Y is out of range, return 0
      #   return(0)
      # }else{
      #   #If Y is in range, find the closest value in the grid returned by density
      #     #This is max value returned by dens such that dens <= y
      #   return(dens[max(which(dens$x <= y))])
      # }

      index <- base::findInterval(y, dens$x)
      y_plus_zeros <- c(0,dens$y,0)
      return(y_plus_zeros[index+1])
    }
    dens1_fn <- function(y){dens_fn(dens1,y)}
    dens2_fn <- function(y){dens_fn(dens2,y)}

    ygrid <- base::seq(from = min(c(dens1$x,dens2$x)),
                       to = max(c(dens1$x,dens2$x)),
                       length.out = 10000)

    #Compute numeric integral of 1/2 |dens1_fn - dens2_fn|
    TV <- 1/2* base::sum( base::abs( dens1_fn(ygrid) - dens2_fn(ygrid) ) * base::diff(ygrid)[1] )
    return(TV)
  }
  else{base::stop("method must be 'density' or 'bins'")}
}



# ##Test TV function
# y1 <- rnorm(10^2, mean = 0, sd = 1)
# y2 <- rnorm(10^2, mean = 1, sd  = 1)
#
# TV_distance_fn(y1,y2)
# TV_distance_fn(y1,y2, method = "bins")
#
# integrate(f = function(.x){ 1/2* abs( dnorm(.x) - dnorm(.x-1) ) } ,
#           lower = -5, upper =5)
