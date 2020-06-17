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

# Used in edm-work/code/simulated/egDeyle/egDeyle.rnw to do movie of how EDM
# works, and to show incorrect rEDM results. Don't want it in package though, so
# maybe just temporarily here, and making plot_explain_edm().
simplexPlot = function(Nx.lags, library, tstar,
                       cobwebbing = TRUE,
                       x.lab = expression("X(t-1)"),
                       y.lab = expression("X(t)"),
                       main = "All the points in lagged space",
                       tstar.col = "blue", tstar.pch = 1, tstar.cex = 1,
                       psivec = NULL, neigh.plot = FALSE, neigh.proj = FALSE,
                       pred.plot = FALSE, pred.rEDM = FALSE, true.val = FALSE,
                       legend.plot = TRUE)
    {
      # Create figure for animation to demonstrate the steps in simplex
      #  projection.
      # Args:
      #   Nx.lags: dataframe with time indexing the row, and columns named
      #         Nt, Ntmin1, Xt, Xtmin1, Xtmin2 and, if E!=0, columns XtPredEeq1,
      #         XtPredEeq2, etc. If the old xt etc. names they get converted
      #         to Xt.
      #   library: dataframe stemming from Nx.lags but with dist from x(tstar)
      #    and resulting rank and weight. Has NA for x(tstar) row.
      #   pred.plot: plot predicted value from manual calculations
      #   pred.rEDM: plot predicted value from rEDM
      #   true.val: plot the true next value
      #   legend.plot: include a legend or not

      # if(tstar.plot)
      #    { toPlot = Nx.lags} else
      #    { toPlot = library}
      Nx.lags = Nx.to.NX.format(Nx.lags) # convert older xt headings to Xt
      par(pty="s")
      plot(Nx.lags$Xtmin1, Nx.lags$Xt,
           xlab = x.lab, ylab = y.lab,
           xlim = Xt.axes.range, ylim = Xt.axes.range,
           main = main)
      points(Nx.lags[c(tstar), c("Xtmin1", "Xt")], col = tstar.col,
              pch = tstar.pch, cex = tstar.cex)
      if(neigh.plot) points(library[psivec, c("Xtmin1", "Xt")],
            pch = 19, col = "red")
      if(neigh.proj)
        {  points(library[psivec+1, c("Xtmin1", "Xt")],
                  pch = 19, col = "purple", cex = 0.5)  # smaller
           for(i in 1:length(psivec))
             {

               iArrows(pull(library[psivec, "Xtmin1"]),
                       pull(library[psivec, "Xt"]),
                       pull(library[psivec+1, "Xtmin1"]),
                       pull(library[psivec+1, "Xt"]),
                       curve = 1, size = 0.7, h.lwd = 2, sh.lwd = 2,
                       width = 1, sh.col = "purple")
               # Example:
               # iArrows(0, 0, 4, 4, curve = 1, size = 0.7, h.lwd = 2,
               #  sh.lwd=2, width = 1, sh.col="red")
             }
        }
      if(pred.plot)
        {
          points(Nx.lags[tstar, "Xt"], XtstarPlus1, col = tstar.col,
            pch = 8)
          iArrows(pull(Nx.lags[tstar, "Xtmin1"]),
                       pull(Nx.lags[tstar, "Xt"]),
                       pull(Nx.lags[tstar, "Xt"]),
                       XtstarPlus1,
                       curve = 1, size = 0.7, h.lwd = 2, sh.lwd = 2,
                       width = 1, sh.col = tstar.col)
          abline(h = pull(library[psivec+1, "Xt"]), col = "lightgrey")
        }
      if(pred.rEDM)
        {
          points(Nx.lags[tstar, "Xt"], Nx.lags[tstar+1, "XtPredEeq2"],
            col = "red", cex = 1.5, lwd = 2)
        }
      if(true.val)
        {
          points(Nx.lags[tstar+1, c("Xtmin1", "Xt")],
                 cex = 1.5, lwd = 2, col = "darkgreen")
        }
      if(legend.plot)
        {
          legend("bottomleft",
                 pch=c(tstar.pch, 19, 8, 1, 1),
                 leg=c("x(t*)", "neighbours", "x(t*+1) pred (wt avge)",
                       "rEDM pred", "true x(t*+1)"),
                 col=c(tstar.col, "red", tstar.col, "red", "darkgreen"),
                 cex=0.85)
        }
   }
