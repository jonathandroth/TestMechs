# test_that("breakdown_defier_share matches compute_min_defiers_binaryM when a binary M is provided", {
#   #Use the fn for the binary M case
#   breakdown_binary <-
#     TestMechs::compute_min_defiers_binaryM(df = TestMechs::burstzyn_data,
#                                            d = "condition2",
#                                            m = "signed_up_number",
#                                            y = "applied_out_fl",
#                                            num_Ybins = 5)$defier_lb
#
#   breakdown_multiple <-
#     TestMechs:::breakdown_defier_share(df = TestMechs::burstzyn_data,
#                                        d = "condition2",
#                                        m = "signed_up_number",
#                                        y = "applied_out_fl",
#                                        num_Ybins = 5)
#
#   #Check that the two are equal
#   expect_equal(breakdown_binary, breakdown_multiple, tolerance = 0.01)
# })

#the function compute_min_defiers_binary has been removed. we should check this more directly by hand
