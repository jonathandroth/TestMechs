#' @title Hypothesis test for the sharp null via Kitagawa (2015)
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   mediator must be discrete.
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param B Bootstrap size, default is 500
#' @param alpha Significance level. Default value is .05
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
test_sharp_null_toru <- function(df, d, m, y, B = 500, alpha = .05, num_Ybins = NULL, cluster = NULL) {

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)


  yvec <- df[[y]]

  num_Ybins <- NULL
  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }

  dvec <- df[[d]]
  mvec <- df[[m]]

  if (is.null(cluster)) {
    clustervec <- NULL
  } else {
    clustervec <- df[[cluster]]
  }
  pval <- mergedZtest(yvec, mvec, dvec, c(0,1),
                      xis = .07,
                      B = 500, alpha = alpha, clustervec = clustervec)$pvals

  return(list(reject = (pval < alpha)))

}

boot_compute_stat_multi <- function(boot_dist, m, lambda, K, L0, L1, limit, xis, xiN){

  temp <- list("P_D0" = boot_dist[1:L0, ], "P_D1" = boot_dist[(L0 + 1):(L0 + L1), ])

  Tw <- compute_stat_multi(temp, m, lambda, K, limit, xis, xiN)

  return(Tw)

}

boot_stats_multi <- function(dists, m, lambda, K, xis, xiN, B, limit, Z_order, data = NULL){

  boot_stats <- matrix(0,xiN,B)

  if (is.null(data[["cluster"]])) {

    L0 <- dim(dists$P_D0)[1]
    L1 <- dim(dists$P_D1)[1]

    YD_dist <- gen_YD_dist(dists, lambda, K)

    for (b in 1:B){
      boot_dist <- gen_boot_dist(YD_dist, m, K, L0, L1)
      boot_stats[, b] <- boot_compute_stat_multi(boot_dist, m, lambda, K, L0, L1, limit, xis, xiN)
    }

  } else {


    for (b in 1:B){
      data_b <- get_bootstrap_draw_clustered(data, "Z", "D", "Y",
                                             cluster = "cluster",
                                             fix_n1 = F)

      m_b <- rep(0, K)
      for (k in 1:K){
        m_b[k] <- nrow(subset(data, data_b$Z == Z_order[k]))
      }
      dists_b <- gen_dists_multi(data_b, m_b, Z_order, K)
      boot_stats[, b] <- compute_stat_multi(dists_b, m_b, lambda, K, limit, xis, xiN)
    }

  }

  return(boot_stats)

}

compute_PQgrids_multi <- function(Q_D0, P_D0, Q_D1, P_D1, limit){

  L1 <- length(P_D1)
  L0 <- length(P_D0)

  ### Constructing intervals for the test statistic with D=1 ###

  if (Q_D1[L1] == 0){
    #Case where there is no mass in the D=1, Z=0 distribution
    p_D1 <- P_D1 - c(0, P_D1[1:(L1 - 1)])
    Q1 <- 0
    P1 <- min(p_D1)
  }
  else{
    D1 <- cbind(P_D1, c(0, P_D1[1:(L1 - 1)]), Q_D1, c(0, Q_D1[1:(L1 - 1)]))
    D1 <- subset(D1, D1[, 3] > D1[ ,4])

    if (dim(D1)[1] > limit){
      #Coarsen the grid in a way that will preserve the maximum
      #Check how many intervals there are where q_D1>p_D1 and p_D1>q_D1
      q_D1 <- D1[, 3] - D1[, 4]
      p_D1 <- D1[, 1] - D1[, 2]
      positive <- sum(q_D1 > p_D1)
      if(Q_D1[L1] < 0.5 & positive > 0 & positive < dim(D1)[1]){
        #Can exclude any point where the difference is negative
        D1 <- subset(D1, q_D1 > p_D1)
      }
      if (positive > 1){
        #Merge consecutive points with p_D1 = 0
        D1 <- cbind(D1, c(0, D1[2:dim(D1)[1], 1] ==  D1[1:(dim(D1)[1] - 1), 2]))
        ind <- D1[, 5] + c(D1[2:dim(D1)[1], 5], 0)

        D1 <- subset(D1, ind < 2)
        D1[1:(dim(D1)[1] - 1), 3] <- (1- D1[2:dim(D1)[1], 5]) * D1[1:(dim(D1)[1] - 1), 3] + D1[2:dim(D1)[1], 5] * D1[2:dim(D1)[1], 3]
        D1 <- subset(D1, D1[, 5] == 0)
      }
    }
    Q1 <- D1[, 3] - matrix(D1[, 4], nrow = dim(D1)[1], ncol = dim(D1)[1], byrow = T)
    Q1 <- Q1[lower.tri(Q1, T)]
    P1 <- D1[, 1] - matrix(D1[, 2], nrow = dim(D1)[1], ncol = dim(D1)[1], byrow = T)
    P1 <- P1[lower.tri(P1, T)]
  }

  ### Constructing intervals for the test statistic with D=0 ###

  if (P_D0[L0] == 0){
    #case where there is no mass in the D=0 Z=1 distribution
    q_D0 <- Q_D0 - c(0, Q_D0[1:(L0 - 1)])
    P0 <- 0
    Q0 <- min(q_D0)
  }
  else{
    D0 <- cbind(P_D0, c(0, P_D0[1:(L0 - 1)]), Q_D0, c(0, Q_D0[1:(L0 - 1)]))
    D0 <- subset(D0, D0[, 1] > D0[, 2])

    if (dim(D0)[1] > limit){
      #Coarsen the grid in a way that will preserve the maximum
      #Check how many intervals there are where q_D0>p_D0 and p_D0>q_D0
      p_D0 <- D0[, 1] - D0[, 2]
      q_D0 <- D0[, 3] - D0[, 4]
      positive <- sum(p_D0 > q_D0)

      if(P_D0[L0] < 0.5 & positive > 0 & positive < dim(D0)[1]){
        #Can exclude any point where the difference is negative
        D0 <- subset(D0, p_D0 > q_D0)
      }

      if (positive > 1){
        #Merge consecutive points with q_D0 = 0
        D0 <- cbind(D0, c(0, D0[2:dim(D0)[1], 3] ==  D0[1:(dim(D0)[1] - 1), 4]))
        ind <- D0[, 5] + c(D0[2:dim(D0)[1], 5], 0)

        D0 <- subset(D0, ind < 2)
        D0[1:(dim(D0)[1] - 1), 1] <- (1- D0[2:dim(D0)[1], 5]) * D0[1:(dim(D0)[1] - 1), 1] + D0[2:dim(D0)[1], 5] * D0[2:dim(D0)[1], 1]
        D0 <- subset(D0, D0[, 5] == 0)
      }
    }
    Q0 <- D0[, 3] - matrix(D0[, 4], nrow = dim(D0)[1], ncol = dim(D0)[1], byrow = T)
    Q0 <- Q0[lower.tri(Q0, T)]
    P0 <- D0[, 1] - matrix(D0[, 2], nrow = dim(D0)[1], ncol = dim(D0)[1], byrow = T)
    P0 <- P0[lower.tri(P0, T)]
  }

  return(list("Q1" = Q1,"P1" = P1, "Q0" = Q0, "P0" = P0))
}

compute_stat_multi <- function(dists, m, lambda, K, limit, xis, xiN){

  Tw <- matrix(0, nrow = K, ncol = xiN)

  for (k in 1:(K - 1)){

    #Force the dists to be matrices (this avoids errors if have 1 row and cast as vector)
    dists$P_D0 <- matrix(dists$P_D0, ncol = K)
    dists$P_D1 <- matrix(dists$P_D1, ncol = K)

    PQ <- compute_PQgrids_multi(dists$P_D0[, k], dists$P_D0[, k + 1] , dists$P_D1[, k], dists$P_D1[, k + 1], limit)

    Q_P_D1 <- PQ$Q1 - PQ$P1
    P_Q_D0 <- PQ$P0 - PQ$Q0

    weight_D1 <- sqrt(lambda[k] * PQ$Q1 * (1 - PQ$Q1) + (1 - lambda[k]) * PQ$P1 * (1 - PQ$P1))
    weight_D0 <- sqrt(lambda[k] * PQ$Q0 * (1 - PQ$Q0) + (1 - lambda[k]) * PQ$P0 * (1 - PQ$P0))

    for (x in 1:xiN){

      xi <- xis[x]

      weight_D1 <- pmax(xi, weight_D1)
      weight_D0 <- pmax(xi, weight_D0)

      weightedQ_P_D1 <- Q_P_D1 / weight_D1
      weightedP_Q_D0 <- P_Q_D0 / weight_D0

      Tw[k, x] <- (as.numeric(m[k]) * as.numeric(m[k + 1]) / (m[k] + m[k + 1])) ^ (0.5) * max(weightedQ_P_D1, weightedP_Q_D0)
    }
  }

  Twfinal <- rep(0,xiN)
  for (x in 1:xiN){
    Twfinal[x] <- max(Tw[, x])
  }

  return(Twfinal)

}

gen_boot_dist <- function(YD_dist, m , K, L0, L1){

  boot_dist <- matrix(0, L0 + L1, K)

  for (k in 1:K){
    boot_data <- sample.int((L0 + L1), size = m[k], replace = TRUE, prob = YD_dist[, k])
    boot_pdf <- tabulate(boot_data, nbins = (L0 + L1))/m[k]
    boot_D0 <- cumsum(boot_pdf[1:L0])
    boot_D1 <- cumsum(boot_pdf[(L0 + 1):(L0 + L1)])
    boot_dist[, k] <- c(boot_D0, boot_D1)
  }

  return(boot_dist)
}

gen_dists_multi <-function(data, m, Z_order, K){

  YZ_D0 <- subset(data, data$D == 0)
  YZ_D0 <- YZ_D0[order(YZ_D0$Y), ]
  YZ_D1 <- subset(data, data$D == 1)
  YZ_D1 <- YZ_D1[order(YZ_D1$Y), ]

  Y0_grid <- unique(YZ_D0$Y)
  Y1_grid <- unique(YZ_D1$Y)

  P_D0 <- matrix(0, nrow = length(Y0_grid), ncol = K)
  P_D1 <- matrix(0, nrow = length(Y1_grid), ncol = K)

  for (k in 1:K){
    Y_D0 <- subset(YZ_D0$Y, YZ_D0$Z == Z_order[k])
    P_D0[, k] <- findInterval(Y0_grid, Y_D0) / m[k]
    Y_D1 <- subset(YZ_D1$Y, YZ_D1$Z == Z_order[k])
    P_D1[, k] <- findInterval(Y1_grid, Y_D1) / m[k]
  }

  return(list("P_D0" = P_D0, "P_D1" = P_D1))
}

gen_YD_dist <- function(dists, lambda, K){

  if(K > 2){
    propscore <- dists$P_D1[dim(dists$P_D1)[1], ]
    propscore_comp <- (propscore[2:(K - 1)]-propscore[1:(K - 2)]) ^ 2 - (propscore[2:(K - 1)]-propscore[3:K]) ^ 2
    propscore_comp <- c(1, propscore_comp > 0, 0)
  }
  else{
    propscore_comp <- c(1, 0)
  }

  dist_sum <- diag(lambda, nrow = K, ncol = (K - 1))
  dist_sum[2:K, ] <- dist_sum[2:K, ] + diag((1 - lambda), nrow = (K - 1), ncol = (K - 1))
  Ydist_D0 <- dists$P_D0 %*% dist_sum
  Ydist_D1 <- dists$P_D1 %*% dist_sum
  Ydist_D0D1 <- rbind(Ydist_D0, Ydist_D1)

  YD_dist <- matrix(0, dim(Ydist_D0D1)[1], K)
  for (k in 1:K){
    YD_dist[, k] <- Ydist_D0D1[, k - 1 + propscore_comp[k]]

  }

  return(YD_dist)
}

mergedZtest <- function(Y, D, Z, Z_order,
                        xis = c(sqrt(0.005 * (1 - 0.005)),
                                sqrt(0.05 * (1 - 0.05)),
                                sqrt(0.1 * (1 - 0.1)),
                                1),
                        B = 500, alpha = c(0.1, 0.05, 0.01),
                        clustervec = NULL){

  data <- data.frame("Y" = Y, "D" = D, "Z" = Z)
  if (!is.null(clustervec)) data[["cluster"]] <- clustervec
  data <- na.omit(data)

  limit <- (10 ^ 4)
  limit <- floor((-1 + sqrt(1 + (8 * limit))) / 2)

  xis <- sort(xis)
  xiN <- length(xis)
  K <- length(Z_order)

  m <- rep(0, K)
  for (k in 1:K){
    m[k] <- nrow(subset(data, data$Z == Z_order[k]))
  }
  lambda <- m[2:K] / (m[1:(K - 1)] + m[2:K])

  dists <- gen_dists_multi(data, m, Z_order, K)

  Tw <- compute_stat_multi(dists, m, lambda, K, limit, xis, xiN)

  boot_stats <- boot_stats_multi(dists, m, lambda, K, xis, xiN, B, limit, Z_order = Z_order,
                                 data = data)

  pval <- rowSums(Tw <= boot_stats) / B
  cvs <- apply(boot_stats, 1, quantile, probs = (1 - alpha))

  return(list("teststat" = Tw, "pvals" = pval, "cvs" = cvs))

}
