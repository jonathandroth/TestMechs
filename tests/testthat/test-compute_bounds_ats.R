# test_that("Lee bounds and ours coincide if don't do correction for point mass", {
#   library(leebounds)
#   leeresults <-
#     leebounds::leebounds(TestMechs::kerwin_data %>%
#                 dplyr::rename(outcome = EL_EGRA_PCA_Index,
#                        selection = primarily_leblango,
#                        treat = treated))
#
#   bounds_ats <-
#   compute_bounds_ats(df = kerwin_data,
#                      d= "treated",
#                      m = "primarily_leblango",
#                      y = "EL_EGRA_PCA_Index", adjust_for_point_mass = F)
#   expect_equal(bounds_ats$lb, leeresults$lower_bound,
#                tolerance = 0.003)
#
#   expect_equal(bounds_ats$ub, leeresults$upper_bound,
#                tolerance = 0.003)
#
# })
#


test_that("Bounds get wider when correcting for point mass", {

  bounds_ats_adj <-
    compute_bounds_ats(df = kerwin_data,
                       d= "treated",
                       m = "primarily_leblango",
                       y = "EL_EGRA_PCA_Index", adjust_for_point_mass = T)
  bounds_ats <-
    compute_bounds_ats(df = kerwin_data,
                       d= "treated",
                       m = "primarily_leblango",
                       y = "EL_EGRA_PCA_Index",
                       adjust_for_point_mass = F)
  expect_lt(bounds_ats_adj$lb-0.002, bounds_ats$lb)

  expect_gt(bounds_ats_adj$ub+0.002, bounds_ats$ub)

})
