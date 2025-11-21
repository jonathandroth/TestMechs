skip_if_not_installed("fixest")
skip_if_not_installed("osqp")
skip_if_not_installed("Rglpk")

library(haven)
library(ggplot2)

setup_baranov_data <- function() {
  data("baranov_data", package = "TestMechs")
  mother_data <- mother_data |> dplyr::filter(THP_sample == 1)
  mother_data_plus_fake <- dplyr::bind_rows(
    mother_data |> dplyr::mutate(fake = 0),
    mother_data |> dplyr::mutate(treat = 0, fake = 1)
  )

  list(
    mother_data = mother_data,
    mother_data_plus_fake = mother_data_plus_fake
  )
}

expect_equivalent_results <- function(ref, ...) {
  lapply(list(...), function(candidate) {
    testthat::expect_equal(candidate, ref, tolerance = 1e-8)
  })
}

expect_distinct_results <- function(lhs, rhs, tolerance = 1e-8) {
  testthat::expect_false(
    isTRUE(all.equal(lhs, rhs, tolerance = tolerance)),
    info = "Specifications with additional controls should differ from no-control estimates"
  )
}

expect_ggplots_equal <- function(..., strict = TRUE) {
  # We only really need svglite here
  testthat::skip_if_not_installed("svglite")

  plots <- list(...)
  if (length(plots) == 1 && is.list(plots[[1]])) {
    plots <- plots[[1]]
  }

  if (length(plots) < 2) {
    stop("expect_ggplots_equal_svg(...) requires at least two ggplot objects.")
  }

  # Helper: render a ggplot to an SVG string
  to_svg <- function(p) {
    if (!inherits(p, "ggplot")) {
      stop("All inputs to expect_ggplots_equal_svg must be ggplot objects (or a list of them).")
    }
    svglite::stringSVG(print(p))
  }

  ref_svg <- to_svg(plots[[1]])

  for (i in 2:length(plots)) {
    svg_i <- to_svg(plots[[i]])

    if (strict) {
      # Strict: SVG strings must be exactly identical
      testthat::expect_equal(
        svg_i,
        ref_svg,
        info = paste("Plot", i, "differs from plot 1")
      )
    } else {
      # Non-strict option if you later want to do something fuzzier
      testthat::expect_true(
        identical(svg_i, ref_svg),
        info = paste("Plot", i, "differs from plot 1")
      )
    }
  }

  invisible(TRUE)
}

test_that("test_sharp_null regression adjustments match baseline (binary M)", {
  datasets <- setup_baranov_data()

  base_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc"
  )

  trivial_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat"
  )

  augmented_res <- test_sharp_null(
    df = datasets$mother_data_plus_fake,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat | fake",
    analytic_variance = TRUE
  )

  expect_equivalent_results(base_res, trivial_res, augmented_res)
})

test_that("test_sharp_null regression adjustments match baseline (nonbinary M)", {
  datasets <- setup_baranov_data()

  base_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc"
  )

  trivial_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat"
  )

  expect_equivalent_results(base_res, trivial_res)
})

test_that("fixest control and IV specifications align for nonbinary M", {
  datasets <- setup_baranov_data()

  controls_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat + factor(interviewer)"
  )

  fe_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat | interviewer"
  )

  iv_res <- test_sharp_null(
    df = datasets$mother_data |> dplyr::mutate(treativ = treat),
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ (treat = treativ) | interviewer"
  )

  expect_equivalent_results(controls_res, fe_res, iv_res)
})

test_that("fixest control and IV specifications align for binary M", {
  datasets <- setup_baranov_data()

  controls_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat + factor(interviewer)"
  )

  fe_res <- test_sharp_null(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ treat | interviewer"
  )

  iv_res <- test_sharp_null(
    df = datasets$mother_data |> dplyr::mutate(treativ = treat),
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    method = "CS",
    num_Ybins = 5,
    cluster = "uc",
    reg_formula = "~ (treat = treativ) | interviewer"
  )

  expect_equivalent_results(controls_res, fe_res, iv_res)
})

test_that("lb_frac_affected respects regression choices", {
  datasets <- setup_baranov_data()

  baseline <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    num_Ybins = 5,
    allow_min_defiers = TRUE
  )

  trivial <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat",
    allow_min_defiers = TRUE
  )

  controls <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat + factor(interviewer)",
    allow_min_defiers = TRUE
  )

  fe <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat | interviewer",
    allow_min_defiers = TRUE
  )

  iv <- lb_frac_affected(
    df = datasets$mother_data |> dplyr::mutate(treativ = treat),
    d = "treat",
    m = "relationship_husb",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ (treat=treativ) | interviewer",
    allow_min_defiers = TRUE
  )

  expect_equivalent_results(baseline, trivial)
  expect_equivalent_results(controls, fe, iv)

  expect_distinct_results(baseline, controls)
  expect_distinct_results(trivial, controls)
})

test_that("lb_frac_affected matches across mediator vector regressions", {
  datasets <- setup_baranov_data()

  baseline <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = c("relationship_husb", "grandmother"),
    y = "motherfinancial",
    num_Ybins = 5,
    allow_min_defiers = TRUE
  )

  trivial <- lb_frac_affected(
    df = datasets$mother_data,
    d = "treat",
    m = c("relationship_husb", "grandmother"),
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat",
    allow_min_defiers = TRUE
  )

  expect_equivalent_results(baseline, trivial)
})


test_that("compute_bounds_ats_new agrees across reg specs", {
  baseline <- compute_bounds_ats_new(df = datasets$mother_data, d = "treat", y = "motherfinancial", m= "relationship_husb", at_group = 5)
  trivial <- compute_bounds_ats_new(df = datasets$mother_data, d = "treat", y = "motherfinancial", m= "relationship_husb", at_group = 5, reg_formula = "~ treat")
  expect_equivalent_results(baseline, trivial)
})


test_that("compute_bounds_ats_new agrees across reg specs", {
  controls <- compute_bounds_ats_new(df = datasets$mother_data, d = "treat", y = "motherfinancial",
                                     m= "relationship_husb", at_group = 5,
                                     reg_formula = "~ treat + factor(interviewer)")
  fes <- compute_bounds_ats_new(df = datasets$mother_data, d = "treat", y = "motherfinancial",
                                m= "relationship_husb", at_group = 5,
                                reg_formula = "~ treat | interviewer")
  iv <- compute_bounds_ats_new(df = datasets$mother_data |> dplyr::mutate(treativ = treat), d = "treat", y = "motherfinancial",
                               m= "relationship_husb", at_group = 5,
                               reg_formula = "~ (treat=treativ) | factor(interviewer)")
  expect_equivalent_results(controls, fes, iv)
})

test_that("partial_density_plot agrees across regression specifications", {
  datasets <- setup_baranov_data()


  plot_controls <- partial_density_plot(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat + factor(interviewer)"
  )

  plot_fe <- partial_density_plot(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat | interviewer"
  )

  expect_ggplots_equal(plot_controls, plot_fe)
})



test_that("partial_density_plot agrees across regression specifications", {
  datasets <- setup_baranov_data()

  plot_base <- partial_density_plot(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    num_Ybins = 5
  )

  plot_trivial <- partial_density_plot(
    df = datasets$mother_data,
    d = "treat",
    m = "grandmother",
    y = "motherfinancial",
    num_Ybins = 5,
    reg_formula = "~ treat"
  )

  expect_ggplots_equal(plot_base, plot_trivial)
})
