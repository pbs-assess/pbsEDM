# Nx_lags_orig3.R - rerunning with latest version of rEDM. Currently had 0.7.4
# saved, now getting 1.7.3 from CRAN, dated 17 Dec 2020. Now redoing with
# Simplex() since simplex() seems deprecated, and saving these as _3 (simplex
# was _2, older ones have no number).

# Generate the example N time series from Carrie's original function (that
#  demonstrates the concerns with rEDM), differences and lagged values, and
#  predictions from rEDM and from my original code.
#  Andy running once (no need for others
#  to run, hence the hardwired pathnames), based on
#  edm-work/code/simulated/egDeyle/tStarLoop/tStarLoop19.rnw.
#  See code at end for updating column names with new notation.

rm(list=ls())
load_all()
# source("../../edm-work/code/functions.r")
source("../../edm-work/code/simulated/sockeye-simulated/SockeyeSim.r")

set.seed(42)
T = 100                   # Number of years of simulated data
simulated = salmonTraj(nyears = T+5)   # Need +5 to get T simulated spawners
                         # Simulated annual spawner abundances and recruitments
                         #  (as a list object)
N = simulated$S          # Just the simulated spawners
tvec = 1:T

# First-difference the data
# Original data are
# N(1), N(2), N(3), ..., N(T),
# First difference:
# X(t) = N(t+1) - N(t)  for t=1, ..., T-1.
X = N[-1] - N[-length(N)]    # X = first-difference vector of
                             #       X(1), X(2), ..., X(T-1)
# expect_equal(X, simple_ts)  is true, using my originally saved simple_ts
# Not standardising since univariate.

# Create Nx.lags with extra columns to represent lagged variables (manually,
# TODO put a test comparing with Luke's function).
# TODO Change Xt to Y_t, etc., check with what gets saved
Nx.lags = data.frame("t" = tvec, "Nt" = N)   # now adding t
Nx.lags = dplyr::tbl_df(Nx.lags)   # TODO change to tibble
Nx.lags = dplyr::mutate(Nx.lags,
                        "Ntmin1" = c(NA, Nt[-T]),
                        "Xt" = c(Nt[-1] - Nt[-T], NA),
                        "Xtmin1" = c(NA, Xt[-T]),
                        "Xtmin2" = c(NA, NA, Xt[-c(T-1, T)])
                        )
Efix = 2

# From tStarLoop19.rnw in Sept 2019:
# Those results (except full prediction values as hard to compare by eye), are
#  the same as for old rEDM, except that rho is now higher, mae is lower, rmse
#  is lower, p-val differenct (but still essentially zero); all const_ values
#  are the same as to be expected. In my 2017 summary I said I calculated rho to
#  be 0.70, so looks like rEDM value now agrees. But see later for exact points...

# Try with new Simplex(), not simplex(), function, need to make a dataframe
#  first:
input = dplyr::select(Nx.lags, t, Xt) %>%
  dplyr::rename(Time = t)

simp.Efix = rEDM::Simplex(dataFrame = input,
                          columns = "Xt",
                          target = "Xt",
                          lib = "1 99",
                          pred = "1 99",
                          E = Efix,
                          verbose = TRUE) #stats_only = FALSE)
stats = rEDM::ComputeError(simp.Efix$Observations,
                           simp.Efix$Predictions) # computes all three stats
rEDM.rho.3 = stats$rho
rEDM.rho.3


# rEDM.points = simp.Efix[,"model_output"][[1]]
rEDM.points = simp.Efix
rEDM.points = tibble::as_tibble(rEDM.points) %>%
  dplyr::rename(obs = Observations,
                pred = Predictions,
                pred_var = Pred_Variance) %>%  # stick with old names
  dplyr::mutate(std_err = sqrt(pred_var),
                pred_min = pred - std_err,
                pred_max = pred + std_err)
                   # min and max but based on std error, not 95% conf intervals
                   # TODO: should change to 2*std_err if using in detail

# mean pred is within the standard errors?
rEDM.points = dplyr::mutate(rEDM.points,
                            in_int = ((obs > pred_min) & (obs < pred_max)))
num.in = sum(rEDM.points$in_int, na.rm=TRUE)       # within interval
num.poss = sum(!is.na(rEDM.points$obs * rEDM.points$pred))  # have obs and pred
percent.in = num.in/num.poss * 100
percent.in
# manually calc absolute error
rEDM.points = dplyr::mutate(rEDM.points,
                            diff = obs - pred)
mae.manual = mean(abs(rEDM.points$diff),
                  na.rm=TRUE)
mae.manual
testthat::expect_equal(mae.manual,
                       stats$MAE)
# That agrees with the rEDM calculated value (as it did for 2017 version).

Nx.lags = dplyr::mutate(Nx.lags,
                        rEDM.pred = c(NA, rEDM.points$pred),
                        rEDM.var = c(NA, rEDM.points$pred_var))
                        # rEDM.pred was XtPredEeq2, but use this to specify it's rEDM
res = EDM_pred_E_2(Nx.lags)   # Andy manual function
Nx.lags = res$Nx.lags
my.full.calcs =  res$my.full.calcs
psi.values = res$psi.values

Nx.lags = dplyr::mutate(Nx.lags,
                        pred.diff = rEDM.pred - my.pred,
                        var.diff = rEDM.var - my.var,
                        pred.ratio = rEDM.pred / my.pred,
                        var.ratio = rEDM.var / my.var)
# round(Nx.lags$pred.diff, 6)
eps = 0.0000001          # how far they have to be apart to investigate
                         # index of the predicted value, t^*+1, is non-zero:
tstarPlus1.big.pred.diff = which(abs(Nx.lags$pred.diff) > eps)
tstarPlus1.big.pred.diff
tstarPlus1.big.var.diff = which(abs(Nx.lags$var.diff) > eps)
tstarPlus1.big.var.diff

# Rename to save them in PBSedm package
# Append with _2 for now to not overwrite original ones - that's using new rEDM
# (Dec 2020) simplex()
# Appending with _3 is using same package but Simplex() function, since that's
# now recommended.
Nx_lags_orig_3 <- Nx.lags
full_calcs_orig_3 =  res$my.full.calcs
psi_orig_3 = res$psi.values

usethis::use_data(Nx_lags_orig_3, overwrite = TRUE)
usethis::use_data(full_calcs_orig_3, overwrite = TRUE)
usethis::use_data(psi_orig_3, overwrite = TRUE)

# These are the ones that are different:
# dplyr::filter(Nx_lags_orig, abs(pred.diff) > 0.0000001)


# This is run after, for the new notation.
NY_lags_example_3 <- Nx_lags_orig_3
names(NY_lags_example_3)[2:6] <- c("N_t", "N_tmin1", "Y_t", "Y_tmin1", "Y_tmin2")
usethis::use_data(NY_lags_example_3, overwrite = TRUE)

# To compare with original:
expect_equal(full_calcs_orig, full_calcs_orig_3)   # my results haven't changed,
                                        # as expected
expect_equal(Nx_lags_orig_2, Nx_lags_orig_3)   # simplex() and Simplex() give
  # different results!
# expect_equal(Nx_lags_orig, Nx_lags_orig_2)   # original and new simplex() match
  # NOT NOW?!?! Do in Nx_lags_orig_2.R, else confusing.
# simplex and Simplex differ:
expect_equal(Nx_lags_orig, Nx_lags_orig_3)
#Error: `Nx_lags_orig` not equal to `Nx_lags_orig_3`.
#Component "rEDM.pred": Mean relative difference: 0.1494412
#Component "rEDM.var": Mean relative difference: 0.4327178
#Component "pred.diff": Mean relative difference: 1
#Component "var.diff": Mean relative difference: 1
#Component "pred.ratio": Mean relative difference: 0.1569088
#Component "var.ratio": Mean relative difference: 0.6525705

expect_equal(NY_lags_example, NY_lags_example_3)
#Error: `NY_lags_example` not equal to `NY_lags_example_3`.
#Component "rEDM.pred": 'is.NA' value mismatch: 2 in current 3 in target
#omponent "rEDM.var": 'is.NA' value mismatch: 2 in current 3 in target
#Component "pred.diff": 'is.NA' value mismatch: 2 in current 3 in target
#Component "var.diff": 'is.NA' value mismatch: 2 in current 3 in target
#Component "pred.ratio": 'is.NA' value mismatch: 2 in current 3 in target
#Component "var.ratio": 'is.NA' value mismatch: 2 in current 3 in target

# For orig versus 2 (I think older rEDM versus new rEDM simplex()) I had
expect_equal(NY_lags_example$rEDM.pred, NY_lags_example_2$rEDM.pred)
#Error: NY_lags_example$rEDM.pred not equal to NY_lags_example_2$rEDM.pred.
#6/100 mismatches (average diff: 0.414)
#[16]  -2.84 - -2.564 == -0.280
#[28]  -2.82 - -2.707 == -0.112
#[40]  -3.82 - -4.389 ==  0.569
#[44]  -2.77 - -2.613 == -0.156
#[82]  -2.55 - -1.591 == -0.954
#[100]   NaN - -0.136 ==    NaN

# But now just get:
expect_equal(NY_lags_example$rEDM.pred, NY_lags_example_3$rEDM.pred)
# [100] NaN - -0.136 == NaN
# and as expected the others will similarly change:
expect_equal(NY_lags_example$rEDM.var, NY_lags_example_3$rEDM.var)
# Error: NY_lags_example$rEDM.var not equal to NY_lags_example_3$rEDM.var.
# 1/100 mismatches
# [100] NaN - 0.00289 == NaN
expect_equal(NY_lags_example$pred.diff, NY_lags_example_3$pred.diff)
# [100] NaN - -2.78e-17 == NaN    # rEDM.pred and my.pred [100] agree (but
# previously the former was NaN), so get these expected results, just confirming.
expect_equal(NY_lags_example$var.diff, NY_lags_example_3$var.diff)
# [100] NaN - 0 == NaN
expect_equal(NY_lags_example$pred.ratio, NY_lags_example_3$pred.ratio)
# [100] NaN - 1 == NaN
expect_equal(NY_lags_example$var.ratio, NY_lags_example_3$var.ratio)
# [23]  3.8e+21 - 3.8e+21 == -524288
# [100]     NaN - 1.0e+00 ==     NaN
# That first one is a numerical thing, don't worry.

expect_equal(psi_orig  , psi_orig_3)

# So my calcs haven't changed, but some of the rEDM ones have for simplex(), but
# not Simplex(), so this not needed now.
# NY_lags_example_2[c(16, 28, 40, 44, 82, 100),]
# Gives the values that have changed, but they don't match mine, though value
# 100 looks to be there now.

#  Now look at  vignette inclusion_issue_3.Rmd
