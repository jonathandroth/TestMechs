test_that("TV returns number close to zero for same distribution", {
  set.seed(0)
  y1 <- rnorm(10^4, mean = 0, sd = 1)
  y2 <- rnorm(10^4, mean = 0, sd  = 1)
  expect_equal(TV_distance_fn(y1,y2),
               0,
               tolerance = 0.03)
})

test_that("TV returns number close to zero for same distribution using bins method", {
  set.seed(0)
  y1 <- rnorm(10^4, mean = 0, sd = 1)
  y2 <- rnorm(10^4, mean = 0, sd  = 1)
  expect_equal(TV_distance_fn(y1,y2,method = "bins", numbins = 40),
               0,
               tolerance = 0.05)
})


test_that("TV returns number close to 1 for very different distributions", {
  set.seed(0)
  y1 <- rnorm(10^4, mean = 0, sd = 1)
  y2 <- rnorm(10^4, mean = 20, sd  = 1)
  expect_equal(TV_distance_fn(y1,y2),
               1,
               tolerance = 0.03)
})

test_that("TV returns number close to 1 for very different distributions using bin method", {
  set.seed(0)
  y1 <- rnorm(10^4, mean = 0, sd = 1)
  y2 <- rnorm(10^4, mean = 20, sd  = 1)
  expect_equal(TV_distance_fn(y1,y2, method = "bins"),
               1,
               tolerance = 0.03)
})

test_that("TV bound close to 1 w/very large treatment effect", {
  expect_equal(compute_tv_ats(df = TestMechs::kerwin_data %>%
                                dplyr::mutate(EL_EGRA_PCA_Index =
                                                ifelse(treated,EL_EGRA_PCA_Index + 20,EL_EGRA_PCA_Index)),
                              d= "treated",
                              m = "primarily_leblango",
                              y = "EL_EGRA_PCA_Index"),
               1,
               tolerance = 0.03)
})


test_that("TV bound close to 0 w/no treatment effect", {
    kerwin_data$fake_y <- rnorm(n = nrow(kerwin_data))
    expect_equal(compute_tv_ats(df = kerwin_data,
                              d= "treated",
                              m = "primarily_leblango",
                              y = "fake_y"),
               0,
               tolerance = 0.03)
})
