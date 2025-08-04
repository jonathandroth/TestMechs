

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
                                         num_gridpoints = 10^5){

  # This is the function F^{-1}(u) := min{y : f(y) >= u}
  inv_cdf <- function(u) {
    min(ecdf_table$y[ecdf_table$cdf >= u])
  }

  #If upper = T, we calculate E[F^{-1} | U \in [1-frac,1]]
  #If upper = F, we calculate E[F^{-1} | U \in [0,frac]]
  #If upper = NULL, we calculate both
  if(is.null(upper)){
    return(
      list(lb =
             trimmed_expectation_from_cdf(
                                          ecdf_table = ecdf_table,
                                          frac = frac,
                                          upper = FALSE,
                                          num_gridpoints = num_gridpoints),
           ub = trimmed_expectation_from_cdf(
                                          ecdf_table = ecdf_table,
                                          frac = frac,
                                          upper = TRUE,
                                          num_gridpoints = num_gridpoints)))
  }

  if(upper == T){
    u_grid <- seq(from = 1-frac, to = 1, length.out = num_gridpoints)
  }else{
    u_grid <- seq(from = 0, to = frac, length.out = num_gridpoints)
  }


  trimmed_expectation <- mean(sapply(u_grid, inv_cdf))

  return(trimmed_expectation)
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
                               max_defier_share = NULL,
                               reg_formula = NULL,
                               num_gridpoints = 10^5){

  #features not yet implemented. Chen, you can fill these in and remove these errors
  if(at_group != 1){stop("function is currently only implemented for at_group = 1")}
  if(!is.null(reg_formula)){stop("reg_formula not yet implemented")}
  if(!is.null(max_defier_share)){stop("function is currently only implemented when max_defier_share = NULL") }

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]


  #compute type shares
  at_share <- base::mean(mvec[dvec == 0])
  nt_share <- base::mean(1-mvec[dvec == 1])
  c_share <- 1 - at_share - nt_share

  #Share of ATs among M=1,D=1
  at_share_m1d1 = at_share / (c_share + at_share)

  #Share of ATs among M=1,D=0
  at_share_m1d0 = 1

  ##Compute bounds on Y(1,1) | G=AT

  #First we compute the CDF of Y | M=1,D=1
  yvalues <- sort(unique(yvec))

  ecdf_y_m1d1_fn <- stats::ecdf(yvec[ mvec ==1 & dvec ==1 ]) #this returns a fn y -> CDF
  ecdf_y_m1d1 <- ecdf_y_m1d1_fn(yvalues) #this is a vector of CDFs evaluated at yvalues


  #Now we compute the bounds on E[Y(1,1) | G=11]
  # The ub is the top at_share_m1d1 of the Y|M=1,D=1 distribution
  # The lb is the bottom at_share_m1d1 of the Y|M=1,D=1 distribution
  bounds_y1 <- trimmed_expectation_from_cdf(ecdf_table = data.frame(y=yvalues, cdf = ecdf_y_m1d1),
                                        frac = at_share_m1d1,
                                        upper = NULL)


  #Similarly we compute the CDF of Y | M=1,D=0
  ## In this case, this trimming argument is trivial, but setting it up this week for generalization
  ecdf_y_m1d0_fn <- stats::ecdf(yvec[ mvec ==1 & dvec ==0 ]) #this returns a fn y -> CDF
  ecdf_y_m1d0 <- ecdf_y_m1d0_fn(yvalues) #this is a vector of CDFs evaluated at yvalues


  #Now we compute the bounds on E[Y(0,1) | G=11]
  # The ub is the top at_share_m1d0 of the Y|M=1,D=1 distribution
  # The lb is the bottom at_share_m1d0 of the Y|M=1,D=1 distribution
  bounds_y0 <- trimmed_expectation_from_cdf(ecdf_table = data.frame(y=yvalues, cdf = ecdf_y_m1d0),
                                        frac = at_share_m1d0,
                                        upper = NULL)



  ##Now compute bounds on ADE for always-takers
  # LB is lb for y1 minus ub for y0
  lb <- bounds_y1$lb - bounds_y0$ub
  ub <- bounds_y1$ub - bounds_y0$lb

  return(data.frame(lb=lb, ub=ub))
}

