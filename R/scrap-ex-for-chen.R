library(devtools)
library(dplyr)
load("data/burstzyn_data.rda")
load_all()

#Bounds on always-takers

TestMechs::compute_bounds_ats(df = burstzyn_data,
                              d = "condition2",
                              m = "signed_up_number",
                              y = "applied_out_fl")


compute_bounds_ats_new(df = burstzyn_data,
                       d = "condition2",
                       m = "signed_up_number",
                       y = "applied_out_fl")





#Bounds on never-takers by flipping the arguments
bounds_flipped <-
TestMechs::compute_bounds_ats(df = burstzyn_data %>% mutate(m2 = 1-signed_up_number,
                                                            d2 = 1-condition2),
                              d = "d2",
                              m = "m2",
                              y = "applied_out_fl")

data.frame(lb = -bounds_flipped$ub,
           ub = -bounds_flipped$lb)


bounds_flipped_new <-
compute_bounds_ats_new(df = burstzyn_data %>% mutate(m2 = 1-signed_up_number,
                                                            d2 = 1-condition2),
                              d = "d2",
                              m = "m2",
                              y = "applied_out_fl")

data.frame(lb = -bounds_flipped_new$ub,
           ub = -bounds_flipped_new$lb)
