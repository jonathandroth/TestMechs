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



