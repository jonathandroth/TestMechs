library(dplyr)
library(haven)
test_that("lb_frac_affected is zero when (Y,M) independent of D", {
  noeffect_df <- rbind(mother_data %>% filter(treat == 1),
                       mother_data %>% filter(treat == 1) %>% mutate(treat = 0))

  expect_equal(lb_frac_affected(df = noeffect_df,
                                d = "treat",
                                m = "grandmother",
                                y = "motherfinancial",
                                num_Ybins = 5),
               0,
               tolerance = 0.001)

  expect_equal(lb_frac_affected(df = noeffect_df,
                                d = "treat",
                                m = "relationship_husb",
                                y = "motherfinancial",
                                num_Ybins = 5),
               0,
               tolerance = 0.001)

  expect_equal(lb_frac_affected(df = noeffect_df,
                                d = "treat",
                                m = c("grandmother", "relationship_husb"),
                                y = "motherfinancial",
                                num_Ybins = 5),
               0,
               tolerance = 0.001)
})



test_that("lb_frac_affected is 1 when Y is colinear with D", {
  bigeffect_df <- rbind(mother_data %>% filter(treat == 1) %>% mutate(motherfinancial = 1),
                        mother_data %>% filter(treat == 1) %>% mutate(treat = 0) %>% mutate(motherfinancial = 2))

  expect_equal(lb_frac_affected(df = bigeffect_df,
                                d = "treat",
                                m = "grandmother",
                                y = "motherfinancial",
                                num_Ybins = 5),
               1,
               tolerance = 0.001)

  expect_equal(lb_frac_affected(df = bigeffect_df,
                                d = "treat",
                                m = "relationship_husb",
                                y = "motherfinancial",
                                num_Ybins = 5),
               1,
               tolerance = 0.001)
})



test_that("lb_frac_affected gets it right for subgroups with big or no effect of D on Y", {
  #Create new_df such that
  # Y is the same in treated and control groups when M=0
  # Y is completely different in treated and controul groups when M=1
  new_df <- rbind(mother_data %>% filter(treat == 1) %>% mutate(motherfinancial = 1) %>%
                          mutate(motherfinancial = 1 ),
                        mother_data %>% filter(treat == 1) %>% mutate(treat = 0) %>%
                          mutate(motherfinancial = case_when(grandmother == 0 ~ 1,
                                                             grandmother == 1 ~ 2) ))

  expect_equal(lb_frac_affected(df = new_df,
                                d = "treat",
                                m = "grandmother",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 0),
               0,
               tolerance = 0.001)

  expect_equal(lb_frac_affected(df = new_df,
                                d = "treat",
                                m = "grandmother",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 1),
               1,
               tolerance = 0.001)
})




test_that("lb_frac_affected matches analytic formula in binary case", {

  new_df <- remove_missing_from_df(df = burstzyn_data,
                                   d = "condition2",
                                   m= "signed_up_number",
                                   y = "applied_out_fl")

  treat_index <- which(new_df$condition2 == 1)
  control_index <- which(new_df$condition2 == 0)

  frac_ats <- mean(new_df$signed_up_number[ control_index ])
  frac_nts <- mean(1 - new_df$signed_up_number[ treat_index ])
  frac_cs <- 1- frac_ats - frac_nts

  diff_partial_densities_1 <-
    pmax(0,
         mean( new_df$applied_out_fl[treat_index] == 1 & new_df$signed_up_number[treat_index] ==1) -
           mean( new_df$applied_out_fl[control_index] == 1 & new_df$signed_up_number[control_index] ==1),
         mean( new_df$applied_out_fl[treat_index] == 0 & new_df$signed_up_number[treat_index] ==1) -
           mean( new_df$applied_out_fl[control_index] == 0 & new_df$signed_up_number[control_index] ==1),
         mean(new_df$signed_up_number[treat_index] ==1) -
           mean(new_df$signed_up_number[control_index] ==1))

  lb_ats <-pmax(0,1/frac_ats * (diff_partial_densities_1 - frac_cs))

  lb_ats_pkg <- lb_frac_affected(df = new_df,
                                 d = "condition2",
                                 m = "signed_up_number",
                                 y = "applied_out_fl",
                                 num_Ybins = 2,
                                 at_group = 1)

  expect_equal(lb_ats_pkg,
               lb_ats,
               tolerance = 0.001)


  diff_partial_densities_0 <-
    pmax(0,
         mean( new_df$applied_out_fl[treat_index] == 1 & new_df$signed_up_number[treat_index] ==0) -
           mean( new_df$applied_out_fl[control_index] == 1 & new_df$signed_up_number[control_index] ==0),
         mean( new_df$applied_out_fl[treat_index] == 0 & new_df$signed_up_number[treat_index] ==0) -
           mean( new_df$applied_out_fl[control_index] == 0 & new_df$signed_up_number[control_index] ==0),
         mean(new_df$signed_up_number[treat_index] ==0) -
           mean(new_df$signed_up_number[control_index] ==0))

  lb_nts <-pmax(0,1/frac_nts * (diff_partial_densities_0))

  lb_nts_pkg <- lb_frac_affected(df = new_df,
                                 d = "condition2",
                                 m = "signed_up_number",
                                 y = "applied_out_fl",
                                 num_Ybins = 2,
                                 at_group = 0)

  expect_equal(lb_nts_pkg,
               lb_nts,
               tolerance = 0.001)

  #check that if at_group = NULL, we match the weighted average
  lb_pooled <- frac_ats / (frac_ats + frac_nts) * lb_ats + frac_nts / (frac_ats + frac_nts) * lb_nts

  lb_pooled_pkg <- lb_frac_affected(df = new_df,
                                    d = "condition2",
                                    m = "signed_up_number",
                                    y = "applied_out_fl",
                                    num_Ybins = 2,
                                    at_group = NULL)

  expect_equal(lb_pooled_pkg,
               lb_pooled,
               tolerance = 0.001)

})




test_that("lb_frac_affected matches Lee bounds lower bound w binary M,Y ", {

  lb_tv <- lb_frac_affected(df = burstzyn_data,
                            d = "condition2",
                            m = "signed_up_number",
                            y = "applied_out_fl",
                            num_Ybins = 2,
                            at_group = 0)

  lb_lee <- bounds_ade_ats(df = burstzyn_data,
                           d = "condition2",
                           m = "signed_up_number",
                           y = "applied_out_fl",
                           at_group = 0)[1] %>% as.numeric()

  expect_equal(lb_tv,
               lb_lee,
               tolerance = 0.001)

})


## test that tv bounds with multi-valued M match the binarized versions at the extremes
test_that("tv bounds with multi-valued M match the binarized versions at the extremes", {

  motherdata_binarized <- mother_data %>% mutate(relationship5 = (relationship_husb == 5),
                                                 relationshipnonzero = (relationship_husb > 1))

  expect_equal(lb_frac_affected(df = motherdata_binarized,
                                d = "treat",
                                m = "relationship5",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 1,
                                max_defiers_share = 0.02 ),
               lb_frac_affected(df = mother_data,
                                d = "treat",
                                m = "relationship_husb",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 5,
                                max_defiers_share = 0.02 ),
               tolerance = 0.001)

  expect_equal(lb_frac_affected(df = motherdata_binarized,
                                d = "treat",
                                m = "relationshipnonzero",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 0,
                                max_defiers_share = 0.02 ),
               lb_frac_affected(df = mother_data,
                                d = "treat",
                                m = "relationship_husb",
                                y = "motherfinancial",
                                num_Ybins = 5,
                                at_group = 0,
                                max_defiers_share = 0.02 ),
               tolerance = 0.001)



})



test_that("with binary M, eps change in defiers share leads to change in lb_frac_affected by eps* (1/thetakk) (lb -1)", {
  mother_data_nonmissing <- remove_missing_from_df(df = mother_data,
                                                   d = "treat",
                                                   m = "grandmother",
                                                   y = "motherfinancial")

  frac_nts <- mean(1- mother_data_nonmissing$grandmother[mother_data_nonmissing$treat == 1])

  eps <- 0.0001

  lb <- lb_frac_affected(df = mother_data,
                         d = "treat",
                         m = "grandmother",
                         y = "motherfinancial",
                         num_Ybins = 5,
                         at_group = 0)

  lb_eps <- lb_frac_affected(df = mother_data,
                             d = "treat",
                             m = "grandmother",
                             y = "motherfinancial",
                             num_Ybins = 5,
                             at_group = 0,
                             max_defiers_share = eps)

  expect_equal( (lb - lb_eps) / eps ,
               1/frac_nts *(1-lb),
               tolerance = 0.001)
})
