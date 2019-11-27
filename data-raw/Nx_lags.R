# Generate the example N time series from Carrie's original function (that
#  demonstrates the concerns with rEDM), differences and lagged values, and
#  predictions from rEDM and from my original code.
#  Andy running once (no need for others
#  to run, hence the hardwired pathnames), based on
#  edm-work/code/simulated/egDeyle/tStarLoop/tStarLoop19.rnw.

rm(list=ls())
source("../../edm-work/code/functions.r")
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
Nx.lags = data.frame("t" = tvec, "Nt" = N)   # now adding t
Nx.lags = dplyr::tbl_df(Nx.lags)   # TODO change to tibble
Nx.lags = dplyr::mutate(Nx.lags,
                        "Ntmin1" = c(NA, Nt[-T]),
                        "Xt" = c(Nt[-1] - Nt[-T], NA),
                        "Xtmin1" = c(NA, Xt[-T]),
                        "Xtmin2" = c(NA, NA, Xt[-c(T-1, T)])
                        )
Efix = 2
simp.Efix = rEDM::simplex(X,
                          E = Efix,
                          stats_only = FALSE)
rEDM.rho = simp.Efix$rho
simp.Efix     # New format for 2019
# From tStarLoop19.rnw in Sept 2019:
# Those results (except full prediction values as hard to compare by eye), are
#  the same as for old rEDM, except that rho is now higher, mae is lower, rmse
#  is lower, p-val differenct (but still essentially zero); all const_ values
#  are the same as to be expected. In my 2017 summary I said I calculated rho to
#  be 0.70, so looks like rEDM value now agrees. But see later for exact points...

rEDM.points = simp.Efix[,"model_output"][[1]]
rEDM.points
rEDM.points = dplyr::tbl_df(rEDM.points)
rEDM.points = dplyr::mutate(rEDM.points,
                            std_err = sqrt(pred_var))
# min and max but based on std error, not 95% conf intervals
rEDM.points = dplyr::mutate(rEDM.points,
                            pred_min = pred - std_err,
                            pred_max = pred + std_err)

# mean pred is within the standard errors?
rEDM.points = dplyr::mutate(rEDM.points,
                            in_int = ((obs > pred_min) & (obs < pred_max)))
num.in = sum(rEDM.points$in_int, na.rm=TRUE)       # within interval
num.poss = sum(!is.na(rEDM.points$obs * rEDM.points$pred))  # have obs and pred
percent.in = num.in/num.poss * 100

# manually calc absolute error
rEDM.points = dplyr::mutate(rEDM.points,
                            diff = obs - pred)
mae.manual = mean(abs(rEDM.points$diff),
                  na.rm=TRUE)
mae.manual

# That agrees with the rEDM calculated value (as it did for 2017 version).

Nx.lags = dplyr::mutate(Nx.lags,
                        rEDM.pred = c(NA, rEDM.points$pred),
                        rEDM.var = c(NA, rEDM.points$pred_var))
                        # rEDM.pred was XtPredEeq2, but use this to specify it's rEDM
Nx_lags <- Nx.lags

usethis::use_data(Nx_lags, overwrite = TRUE)
