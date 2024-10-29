## # Unit tests for test_sharp_null function to see if it runs under reasonable configurations of parameters

## # install.packages("doParallel")
## # install.packages("devtools")
## # devtools::install_github("soonwookwon/lpinfer")
## # devtools::install_github("asheshrambachan/HonestDiD")
## # devtools::install_github("haluong89-bcn/leebounds")
## # devtools::install_github("jonathandroth/MedBounds")
## library(testthat)
## library(parallel)
## library(foreach)
## library(doParallel)
## library(dplyr)
## library(MedBounds)


## test_that("If test_sharp_null runs under reasonable configurations of the parameters 
##           using our application (for each possible type of test).", {
            
##   testdf = MedBounds::kerwin_data;
  
##   methods = c("ARP", "CS", "FSST");
##   cluster = c("Age"); # avoid binary cluster variable for testing
##   weightmats = c("diag", "avar", "identity"); # Weight matrix used to implement FSST.
##   numbins = 5;
##   rearranges = c(T, F);
##   fix_n1s = c(T, F);
##   lambdas = c("dd", "ndd"); # FSST only
##   use_ncs = c(T, F);
##   analytical_vars = c(T, F); # ARP or CS method
##   defiers_shares = c(0, 0.1); # problem with FSST if not 0
  
##   param_values <- list(
##     Method = methods,
##     Cluster = cluster,
##     Num_Ybins = numbins,
##     Rearrange = rearranges,
##     Fix_n1 = fix_n1s,
##     Lambda = lambdas,
##     Use_nc = use_ncs,
##     Defiers_share = defiers_shares,
##     Weight.matrix = weightmats,
##     Analytic_variance = analytical_vars
##   )
  
##   # Parallel Computing setup
  
##   param_df <- expand.grid(param_values)
  
##   temppp <- param_df |> 
##     filter(!(Method == "FSST" & Weight.matrix == "avar"))
##   # # Detect the number of cores
##   # num_cores <- detectCores()
##   # 
##   # # Create a cluster
##   # cl <- makeCluster(num_cores - 1)
##   # 
##   # registerDoParallel(cl)
##   # 
##   # start <- Sys.time()
  
##   # Running the tests
##   for(i in 7:nrow(temppp)) {
##     params <- temppp[i, ]
    
##     if (params$Method == "FSST") {
      
##       if (params$Fix_n1 == FALSE){
        
##         test_sharp_null_result <-
##           MedBounds::test_sharp_null(df = testdf,
##                                      d = "treated",
##                                      m = "primarily_leblango",
##                                      y = "EL_EGRA_PCA_Index",
##                                      method = params$Method,
##                                      cluster = params$Cluster,
##                                      num_Ybins = params$Num_Ybins,
##                                      rearrange = params$Rearrange,
##                                      fix_n1 = FALSE,
##                                      use_nc = params$Use_nc,
##                                      # max_defiers_share = params$Defiers_share, # problem to be fixed
##                                      weight.matrix = params$Weight.matrix,
##                                      lambda = params$Lambda)
        
##         testthat::expect_no_error(test_sharp_null_result)
        
##       }else{ # if Fix_n1 == TRUE, no cluster argument, unless the cluster variable is consistent with treatment (either all or none treated)
        
##         test_sharp_null_result <-
##           MedBounds::test_sharp_null(df = testdf,
##                                      d = "treated",
##                                      m = "primarily_leblango",
##                                      y = "EL_EGRA_PCA_Index",
##                                      method = params$Method,
##                                      num_Ybins = params$Num_Ybins,
##                                      rearrange = params$Rearrange,
##                                      fix_n1 = TRUE,
##                                      use_nc = params$Use_nc,
##                                      # max_defiers_share = params$Defiers_share, # problem to be fixed
##                                      weight.matrix = params$Weight.matrix,
##                                      lambda = params$Lambda)
        
##         testthat::expect_no_error(test_sharp_null_result)
##       }
      
##     } else { # methods CS or ARP
      
##       if(params$Fix_n1 == FALSE){
        
##         test_sharp_null_result <-
##           MedBounds::test_sharp_null(df = testdf,
##                                      d = "treated",
##                                      m = "primarily_leblango",
##                                      y = "EL_EGRA_PCA_Index",
##                                      method = params$Method,
##                                      cluster = params$Cluster,
##                                      num_Ybins = params$Num_Ybins,
##                                      rearrange = params$Rearrange,
##                                      fix_n1 = FALSE,
##                                      use_nc = params$Use_nc,
##                                      max_defiers_share = params$Defiers_share,
##                                      analytic_variance = params$Analytic_variance)
        
##         testthat::expect_no_error(test_sharp_null_result)
        
##       }else{ # no cluster argument if fix_n1 == TRUE
        
##         test_sharp_null_result <-
##           MedBounds::test_sharp_null(df = testdf,
##                                      d = "treated",
##                                      m = "primarily_leblango",
##                                      y = "EL_EGRA_PCA_Index",
##                                      method = params$Method,
##                                      num_Ybins = params$Num_Ybins,
##                                      rearrange = params$Rearrange,
##                                      fix_n1 = TRUE,
##                                      use_nc = params$Use_nc,
##                                      max_defiers_share = params$Defiers_share,
##                                      analytic_variance = params$Analytic_variance)
        
##         testthat::expect_no_error(test_sharp_null_result)
##       }
      
##     }
##   }
  
##   # Stop the cluster
##   # stopCluster(cl)
  
##   print(Sys.time() - start)
  
## # 1. Error in { : task 39 failed - "ℹ In index: 1.
## # Caused by error in `dplyr::left_join()`:
## # ! `by` must be a (named) character vector, list, `join_by()` result, or NULL, not a <factor> object."

##   # solution: cluster <- as.character(cluster) line 57 of nonparametric bootstrap
  
  
##   # task 99
##   # 2. Error in solve.default(beta.sigma) :                                                                  
##   #   system is computationally singular: reciprocal condition number = 4.89739e-17
  
##   MedBounds::test_sharp_null(df = testdf,
##                              d = "treated",
##                              m = "primarily_leblango",
##                              y = "EL_EGRA_PCA_Index",
##                              method = "FSST",
##                              num_Ybins = 5,
##                              rearrange = T, 
##                              fix_n1 = T, 
##                              use_nc = T,
##                              weight.matrix = "avar", #### problem.....
##                              lambda = "ndd")
  
## })


