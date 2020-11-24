## Descriptions of saved data objects.

##' Simple example time series
##'
##' Example simulation results originally from Carrie's salmonTraj() function,
##' that end up some discrepancies in results between `rEDM` and our own EDM
##' code. Contains first difference values, $X(t) = N(t+1) - N(t)$ where $N(t)$
##' are population numbers from the from a simple Ricker-like model. `X` is not
##' standardised (which is okay as univariate time series).
##' @format Vector of length 99 (since original `N(t)` was for 100, with each
##' element representing $X(t) = N(t+1) - N(t)$.
##'
##' @source Generated from running `data-raw/simple_ts.R`.
"simple_ts"

##' Summary results from a simple simulation model analysed using `rEDM` and Andy's
##' original independent code - now using NY_lags_example with new notation
##'
##' Contains the original times series of populations `Nt`, simple lagged values
##' and differences, and then calculations from `rEDM::simplex()` and Andy's
##' original independent code. Indpendent code originally done to understand
##' EDM, but revealed discrepancies between those results and those of `rEDM()`.
##' See `full_calcs_orig` and `psi_orig` for full results.
##'
##' @format A 100x14 tibble where each row is a time step (year) and columns
##'   are:
##'   * `t`: time step, 1, 2, 3, ..., 100.
##'   * `Nt`: original time series of number of spawners, $N(t)$ at time $t$ from
##'   Carrie's `salmonTraj()` function.
##'   * `Ntmin1`: $N(t-1)$, the number at the previous time step (value for
##'      `t=1` is therefore `NA`).
##'   * `Xt`: First difference $X(t) = N(t+1) - N(t)$; (value for `t=100` is
##'   therefore `NA`).
##'   * `Xtmin1`: $X(t-1)$
##'   * `Xtmin2`: $X(t-2)$
##'   * `rEDM.pred`: predicted value from `rEDM::simplex()` with embedding
##'   dimension $E=2$.
##'   * `rEDM.var`: estimated variance of `rEDM.pred` (from `rEDM::simplex()`).
##'   * `my.pred`: predicted value from Andy's `EDM_pred_E_2()` function (saved
##'   in `PBSedm`)
##'   * `my.var`: estimated variance of `my.pred`
##'   * `pred.diff`: the difference between predictions, i.e. `rEDM.pred - my.pred`
##'   * `var.diff`: the difference between estimated variances, i.e. `rEDM.var - my.var`
##'   * `pred.ratio`: ratio of the two estimated predictions, as `rEDM.pred / my.pred`
##'   * `var.ratio`: ratio of the two estimated variances, as `rEDM.var / my.var`
##'
##' @source Generated from Andy running `data-raw/Nx_lags.R`.
"Nx_lags_orig"

##' Summary results from a simple simulation model analysed using `rEDM` and Andy's
##' original independent code
##'
##' Contains the original times series of populations `N_t`, simple lagged values
##' and differences, and then calculations from `rEDM::simplex()` and Andy's
##' original independent code. This was independent code originally done to understand
##' EDM, but revealed discrepancies between those results and those of `rEDM()`.
##' See `full_calcs_orig` and `psi_orig` for full results.
##'
##' @format A 100x14 tibble where each row is a time step (year) and columns
##'   are:
##'   * `t`: time step, 1, 2, 3, ..., 100.
##'   * `N_t`: original time series of number of spawners, $N_t$ at time $t$ from
##'   Carrie's `salmonTraj()` function.
##'   * `N_tmin1`: $N_{t-1}$, the number at the previous time step (value for
##'      `t=1` is therefore `NA`).
##'   * `Y_t`: First difference $Y_t = N_{t+1} - N_t$; (value for `t=100` is
##'   therefore `NA`).
##'   * `Ytmin1`: $Y_{t-1}$
##'   * `Ytmin2`: $Y_{t-2}$
##'   * `rEDM.pred`: predicted value from `rEDM::simplex()` with embedding
##'   dimension $E=2$.
##'   * `rEDM.var`: estimated variance of `rEDM.pred` (from `rEDM::simplex()`).
##'   * `my.pred`: predicted value from Andy's `EDM_pred_E_2()` function (saved
##'   in `PBSedm`)
##'   * `my.var`: estimated variance of `my.pred`
##'   * `pred.diff`: the difference between predictions, i.e. `rEDM.pred - my.pred`
##'   * `var.diff`: the difference between estimated variances, i.e. `rEDM.var - my.var`
##'   * `pred.ratio`: ratio of the two estimated predictions, as `rEDM.pred / my.pred`
##'   * `var.ratio`: ratio of the two estimated variances, as `rEDM.var / my.var`
##'
##' @source Generated from Andy running `data-raw/Nx_lags.R`.
"NY_lags_example"


##' Full results from a simple simulation model analysed using Andy's
##' original independent code (still has `Xt` notation not `Y_t`)
##'
##' Contains a list of dataframes from running `EDM_pred_E_2` on the simple time series.
##'
##' @format List of 100 tibbles, each component `[[tstar]]` for each
##'   `tstar`. Each tibble is 100x5, and contains the original `Xt` and
##'   `Xtmin1`, and then the `d`, `rank` and `w` for that
##'   `tstar`. Tibble `[[1]]` is NULL since calculations can't be done for
##'   $t^*=1$, e.g. $N(t-1)$ does not exist.
##'
##' @source Generated from Andy running `data-raw/Nx_lags.R`.
"full_calcs_orig"

##' Psi values from a simple simulation model analysed using Andy's
##' original independent code
##'
##' For each \eqn{t^*}, hello contains values of \eqn{\psi_1, \psi_2} and \eqn{\psi_3} which are
##' the time indices of the nearest, second-nearest and third-nearest neighbours
##' in \eqn{\bf{x}}-space (vector \eqn{\bf{x}}, so in lagged difference space). There
##' are three nearest neighbours since embedding dimension \eqn{E=2}.
##'
##' @format Dataframe of size 100x4 with one row for each `tstar`, and columns:
##'   * `tstar`
##'   * `psi1`, `psi2`, `psi3`: values of \eqn{psi_1, \psi_2} and \eqn{\psi_3} for that
##'   \eqn{t^*}.
##'
##' @source Generated from Andy running `data-raw/Nx_lags.R`.
"psi_orig"

##' Predicted and Observed Values from rEDM output
##'
##' @format Dataframe of size 100x4 with one row for each `time`, and columns:
##'   * `time`
##'   * `obs`
##'   * `pred`
##'   * `pred_var`
##'
"rEDM_points"
