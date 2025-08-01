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
                                                                       continuous_Y = continuous_Y)


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
