#Write unit tests for the lb_frac_affected function

#lb_frac_affected should match the output of compute_tv_ats when a single binary M is provided

test_that("lb_frac_affected matches compute_tv_ats when a binary M is provided", {
  #Use the fn for the binary M case
  tv_binary <-
  TestMechs::compute_tv_ats(df = TestMechs::kerwin_data,
                              d= "treated",
                              m = "primarily_leblango",
                              y = "EL_EGRA_PCA_Index")

  tv_multiple <-
    TestMechs::lb_frac_affected(df = TestMechs::kerwin_data,
                              d= "treated",
                              m = "primarily_leblango",
                              y = "EL_EGRA_PCA_Index",
                              at_group = c(TRUE,TRUE))

  #Check that the two are equal
  expect_equal(tv_binary, tv_multiple, tolerance = 0.01)
})


#lb_frac_affected should agree with the analytic formulas when M is an ordered scalar
test_that("lb_frac_affected should agree with the analytic formulas when M is an ordered scalar",{
  #Create an M with three values
  df <- TestMechs::kerwin_data
  df$m <- dplyr::ntile(df$primarily_leblango,3)

  #Calculate TV_{22} using analytic formulas
  p_m2_1 <- mean(df$m[df$treated == 1] == 2 )
  p_mgte2_1 <- mean(df$m[df$treated == 1] >= 2 )
  p_mgte2_0 <- mean(df$m[df$treated == 0] >= 2 )

  theta_kk <- p_m2_1 - (p_mgte2_1 - p_mgte2_0)

  max_p_diffs_list <- compute_max_p_difference(dvec = df$treated,
                                               mdf = df[,"m"],
                                               yvec = df$EL_EGRA_PCA_Index,
                                               wvec = rep(1/NROW(df),NROW(df)))

  max_p_diff <- max_p_diffs_list$max_p_diffs[max_p_diffs_list$mvalues == 2]

  analytical_tv2 <- (max_p_diff - (p_mgte2_1-p_mgte2_0))/ theta_kk

  tv_multiple <-
    TestMechs::lb_frac_affected(df = df,
                                         d= "treated",
                                         m = "m",
                                         y = "EL_EGRA_PCA_Index",
                                         at_group = 2)

  expect_equal(analytical_tv2, tv_multiple, tolerance = 0.01)

})


#Test the fractional linear programming wrapper that we wrote
test_that("Fractional linear programming works on simple examples", {

  frac_lp_result <-
  TestMechs:::Rglpk_solve_fractional_LP(obj_numerator = c(1,0),
                            obj_denominator = c(0,1),
                            constant_numerator = 0,
                            constant_denominator = 0,
                            mat = matrix(c(1,0,0,-1),byrow=T,nrow=2),
                            rhs = c(3,-1/2),
                            dir=c("<=","<="),
                            max = TRUE, bounds = NULL
  )

  expect_equal(frac_lp_result$optimum,6, tolerance = 0.01)
  expect_equal(max(abs(frac_lp_result$solution - c(3,1/2))),0, tolerance = 0.01)

  frac_lp_result <-
    TestMechs:::Rglpk_solve_fractional_LP(obj_numerator = c(1,0),
                                          obj_denominator = c(0,1),
                                          constant_numerator = 0,
                                          constant_denominator = 0,
                                          mat = diag(2),
                                          rhs = c(3,1),
                                          dir=c("<=",">="),
                                          max = TRUE, bounds = NULL
    )

  expect_equal(frac_lp_result$optimum, 3, tolerance = 0.01)

  frac_lp_result <- TestMechs:::Rglpk_solve_fractional_LP(obj_numerator = c(-2,1),
                                                          obj_denominator = c(1,3),
                                                          constant_numerator = 2,
                                                          constant_denominator = 4,
                                                          mat = matrix(c(-1,1,2,1),byrow=T,nrow=2),
                                                          rhs = c(4,14),
                                                          dir=c("<=","<="),
                                                          bounds = list(lower=list(ind = c(1,2), val = c(0,0)),
                                                                        upper=list(ind = c(1,2), val = c(Inf,6))),
                                                          max = F)
  expect_equal(frac_lp_result$optimum, -12/11, tolerance = 0.01)
  expect_equal( max(abs(frac_lp_result$solution- c(7,0))), 0, tolerance = 0.01)
})

test_that("get error when monotonicity is violated",{
  expect_error(lb_frac_affected(
    kerwin_data %>% dplyr::mutate(minus_treated = 1-treated),
    d = "minus_treated",
    m = "primarily_leblango",
    y = "EL_EGRA_PCA_Index"))
})



test_that("breakdown_defier_share matches compute_min_defiers_binaryM when a binary M is provided", {
  #Use the fn for the binary M case
  breakdown_binary <-
    TestMechs::compute_min_defiers_binaryM(df = TestMechs::burstzyn_data,
                                           d = "condition2",
                                           m = "signed_up_number",
                                           y = "applied_out_fl",
                                           num_Ybins = 5)$defier_lb

  breakdown_multiple <-
    TestMechs:::breakdown_defier_share(df = TestMechs::burstzyn_data,
                                           d = "condition2",
                                           m = "signed_up_number",
                                           y = "applied_out_fl",
                                           num_Ybins = 5)

  #Check that the two are equal
  expect_equal(breakdown_binary, breakdown_multiple, tolerance = 0.01)
})
