# Some earlier functions that may be used for testing, but will be superceded by
# more general functions. EDM_pred_E_2() was originally in edm-work/code as
# myEDMpred().
# New comment
##' Simplex prediction for embedding dimension E=2
##'
##' Simple code to estimate X(t^*+1) and of its variance, for
##'  all valid values of t^* for a tbl_df Nx.lags. Only for embedding dimension
##'  E=2.
##'
##' @param Nx.lags tbl_df with time indexing the row, and columns
##'   Nt, Ntmin1, Xt, Xtmin1, Xtmin2, rEDM.pred and rEDM.var (though only Xt and
##'   Xtmin1 are used for E=2).
##' @param Efix Embedding dimension, currently only works for `Efix`=2.
##' @return List containing
##'  * Nx.lags: dataframe that was inputted but now now also has columns
##'   `my.pred` and `my.var` for my manual calculations of the predicted value
##'   of $X(t^*+1)$ obtained by omitting $X(t^*)$, and its variance.
##'  * my.full.calcs: list of dataframes, one component `[[tstar]]` for each
##'   `tstar`, that contains full `Xt`, `Xtmin1`, `d`, `rank` and `w` for that tstar
##'  * psi.values: dataframe of values of `psi.vec` (so `psi_1`, `psi_2`, and
##'   `psi_3`, since `Efix`=2), one row for each `t*`.
##' @export
##' @author Andrew Edwards
EDM_pred_E_2 <-  function(Nx.lags,
                          Efix = 2)
  {
  lib.orig = dplyr::select(Nx.lags, Xt, Xtmin1)  # Each row is time t
  usable = lib.orig$Xt * lib.orig$Xtmin1  # to give NA if either unavailable
  tstar.valid = which(!is.na(usable))     # valid values of tstar (not first or
                                          #  last for E=2)
  my.pred = rep(NA, nrow(Nx.lags))        # predicted values to get filled in
  my.var  = rep(NA, nrow(Nx.lags))        #  in for loop:
  my.full.calcs = list()                  # all neighbours, weights, etc. for
                                          #  each value of t^*
  psi.values = data.frame("tstar" = 1:nrow(Nx.lags),
                          "psi1"=NA, "psi2"=NA, "psi3"=NA)
                                          # the indices of the nearest neighbours
                                          #  for each value of tstar
  for(tstar.ind in tstar.valid)
    {
      lib.temp = lib.orig
      lib.temp[tstar.ind+1, "Xt"] = NA  # This is the value we don't know
      lib.temp[tstar.ind+2, "Xtmin1"] = NA

      X.tstar.ind = as.numeric(dplyr::select(lib.temp[tstar.ind,], Xt, Xtmin1))
      # Calculate distance from {\bf x}(t^*) = (X(t^*), X(t^*-1))
      lib.temp = dplyr::mutate(lib.temp,
          d=sqrt((Xt - X.tstar.ind[1])^2 + (Xtmin1 - X.tstar.ind[2])^2))
      # Any that have d = NA cannot then be the result of a projection:
      dIsNA = which(is.na(lib.temp$d))
      dIsNA = dIsNA[dIsNA > 1]       # else get 0 row in next line;
      lib.temp[dIsNA-1, "d"] = NA
      # Rank the valid ones with repsect to distance from x(t^*)
      lib.temp = dplyr::mutate(lib.temp, rank = dplyr::dense_rank(d))

      psivec = rep(NA, Efix+1)
      for(i in 1:length(psivec))
        {
          psivec[i] = which(lib.temp$rank == i)
          psi.values[tstar.ind, paste0("psi", i)] = psivec[i]
        }
      dxstarpsi1 = as.numeric(lib.temp[psivec[1], "d"])  # d of closest neighbour
               # weights for all points though only need Efix+1 closest:
      lib.temp = dplyr::mutate(lib.temp, w = exp(- d / dxstarpsi1))
      sumw = sum(lib.temp[psivec, "w"])
      Xtstar.ind.Plus1 = sum( lib.temp[psivec, "w"] *
                                 lib.temp[psivec+1, "Xt"] ) / sumw
      # Variance:
      Xtstar.ind.Plus1.var = sum( lib.temp[psivec, "w"] *
                     ( lib.temp[psivec+1, "Xt"] - Xtstar.ind.Plus1 )^2 ) / sumw
      my.pred[tstar.ind+1] = Xtstar.ind.Plus1   # the prediction of X(t^*+1)
      my.var[tstar.ind+1] = Xtstar.ind.Plus1.var
      my.full.calcs[[tstar.ind]] = lib.temp  # the nearest neighbours etc. for
                                             #  t^*
    }
  Nx.lags = dplyr::mutate(Nx.lags, my.pred = my.pred, my.var = my.var)
  return(list("Nx.lags" = Nx.lags, "my.full.calcs" = my.full.calcs,
              "psi.values" = psi.values))
  }
