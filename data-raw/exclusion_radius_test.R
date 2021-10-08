# exclusion_radius_test.R - test alternative exclusion radii in rEDM, for use in
#  inclusion issue vignette. Saves the results in pbsEDM (to avoid it requiring rEDM).
#  Run this manually (line-by-line) since will create errors from testthat -
#   check they agree (especially when updating rEDM).

rm(list=ls())
load_all()

input_obs <- dplyr::select(NY_lags_example_3, t, Y_t) %>%
  dplyr::rename(Time = t)

check_rEDM_excl_1 <- rEDM::Simplex(dataFrame = input_obs,
                                   columns = "Y_t",
                                   target = "Y_t",
                                   lib = "1 99",
                                   pred = "1 99",
                                   E = 2,
                                   verbose = TRUE,
                                   exclusionRadius = 1)$Predictions
check_rEDM_excl_1 <- c(NA,
                       check_rEDM_excl_1)   # Needs extra NA to get indexing correct

testthat::expect_equal(NY_lags_example_3$rEDM.pred,
                       check_rEDM_excl_1)
# Fails:
# 2/100 mismatches (average diff: 0.969)
# [76]  1.368 - 0.838 ==  0.53
# [77] -0.381 - 1.027 == -1.41

# So two differences using Simplex() with exclusionRadius=0 (default) and 1.

# From inclusion issue vignette, need:
pbs_calc <- pbsEDM(NY_lags_example_3,
                   lags = list(Y_t = c(0:1))) # A tibble (contains lists)
pbs_pred <- pbs_calc$X_forecast[-length(pbs_calc$X_forecast)]


testthat::expect_equal(pbs_pred,
                       check_rEDM_excl_1)
# Fails:
# 6/100 mismatches (average diff: 0.274)
# [19]  1.0614 - 1.0614 ==  5.82e-06 #small
# [23]  1.4783 - 1.4783 ==  8.08e-06 #small
# [77] -0.3812 - 1.0266 == -1.41e+00
# [90]  2.1994 - 2.1994 ==  1.43e-05 #small
# [95]  0.4123 - 0.1774 ==  2.35e-01
# [98]  0.0799 - 0.0799 ==  3.28e-06 #small
# So [76] now matches pbs_pred, [77] for some reason has changed, [90] is
# probably close enough, but [95] still not the same as pbs_pred.

usethis::use_data(check_rEDM_excl_1,
                  overwrite = TRUE)

# Now try exclusion radius of 2:

check_rEDM_excl_2 <- rEDM::Simplex(dataFrame = input_obs,
                                   columns = "Y_t",
                                   target = "Y_t",
                                   lib = "1 99",
                                   pred = "1 99",
                                   E = 2,
                                   verbose = TRUE,
                                   exclusionRadius = 2)$Predictions
check_rEDM_excl_2 <- c(NA,
                       check_rEDM_excl_2)   # Needs extra NA to get indexing correct

testthat::expect_equal(NY_lags_example_3$rEDM.pred,
                       check_rEDM_excl_2)
# 4/100 mismatches (average diff: 0.931)
# [76]  1.368 - 0.838 ==  0.530  # same as exclR=1
# [77] -0.381 - 1.027 == -1.408  # same as exclR=1
# [95]  0.177 - 0.412 == -0.235
# [97] -0.798 - 0.755 == -1.552
# So this has two new differences from 0

testthat::expect_equal(pbs_pred,
                       check_rEDM_excl_2)
# [19]  1.0614 - 1.0614 ==  5.82e-06 # small
# [23]  1.4783 - 1.4783 ==  8.08e-06 # small
# [77] -0.3812 - 1.0266 == -1.41e+00
# [90]  2.1994 - 2.1994 ==  1.43e-05 # small
# [97] -0.7977 - 0.7547 == -1.55e+00
# [98]  0.0799 - 0.0799 ==  3.28e-06 # small

usethis::use_data(check_rEDM_excl_2,
                  overwrite = TRUE)

# Now try exclusion radius of 3:

check_rEDM_excl_3 <- rEDM::Simplex(dataFrame = input_obs,
                                   columns = "Y_t",
                                   target = "Y_t",
                                   lib = "1 99",
                                   pred = "1 99",
                                   E = 2,
                                   verbose = TRUE,
                                   exclusionRadius = 3)$Predictions
check_rEDM_excl_3 <- c(NA,
                       check_rEDM_excl_3)   # Needs extra NA to get indexing correct

testthat::expect_equal(NY_lags_example_3$rEDM.pred,
                       check_rEDM_excl_3)
# 7/100 mismatches (average diff: 0.64)
# [61]  0.797 - 0.891 == -0.0943
# [64]  0.171 - 0.150 ==  0.0209
# [65]  2.883 - 3.521 == -0.6378
# [76]  1.368 - 0.838 ==  0.5299
# [77] -0.381 - 1.027 == -1.4078
# [95]  0.177 - 0.412 == -0.2350
# [97] -0.798 - 0.755 == -1.5523
# So first three values are different (but weren't for exclR = 2), last four
# have the same differences as for exclR = 2
# (expect to get more differences with larger radius).

testthat::expect_equal(pbs_pred,
                       check_rEDM_excl_3)
# 9/100 mismatches (average diff: 0.413)
# [19]  1.0614 - 1.0614 ==  5.82e-06
# [23]  1.4783 - 1.4783 ==  8.08e-06
# [61]  0.7967 - 0.8910 == -9.43e-02 # new
# [64]  0.1709 - 0.1500 ==  2.09e-02 # new
# [65]  2.8828 - 3.5206 == -6.38e-01 # new
# [77] -0.3812 - 1.0266 == -1.41e+00
# [90]  2.1994 - 2.1994 ==  1.43e-05
# [97] -0.7977 - 0.7547 == -1.55e+00
# [98]  0.0799 - 0.0799 ==  3.28e-06

usethis::use_data(check_rEDM_excl_3,
                  overwrite = TRUE)
