# Some earlier functions that may be used for testing, but will be superceded by
# more general functions. EDM_pred_E_2() was originally in edm-work/code as
# myEDMpred().

# Mostly adapted now in plotting.R, now taken out of analyse_simple... vignette.
##' Plot a movie of data (and EDM results if available) with panels showing data
##'  in various ways
##'
##' Create figure with six panels as a movie. Panels are:
##' N_t vs time and x_t vs time (optional)
##' N_t vs N_{t-1}, x_t vs x_{t-1} (like plotLag2d),
##' x_t vs x_{t-1} vs x_{t-2} (like plotLag3d) and
##' pred vs obs (for given values of E, if provided).
##' All code is in here, may want to parse out into individual plots to be able
##' to more easily reuse.
##' With colours to show the latest points, and star to highlight latest point.
##' Defaults are for `Nx_lags_orig` example.
##'
##' TODO: change options ot have NULl if depend on earlier things. end may need
##' to be one less for default.
##' TODO: change . to _ everywhere
##' TODO: options to save as pdf (maybe - used to be for Sweave but now
##' superceded) or gif to be in a vignette.
##' TODO: adapt to work with Luke's output formats, or just write a wrapper
##' function. Most of this is just looking at data in various ways.
##'
##' @param Nx.lags Dataframe (Tibble???) with time indexing the row, and columns
##'     named Nt, Ntmin1, Xt, Xtmin1, Xtmin2 and, if E!=0, columns XtPredEeq1,
##'     XtPredEeq2, etc. TODO - not same now, using saved example but also want
##'   option for doing Eeq2 etc. from rEDM and from pbsEDM. TODO TODO
##'     This allows Xt to equal Nt if want to plot figures
##'     without the first differencing (basically the 3d plot and the pred
##'     vs obs would be the only new ones). Nt is the original data,
##'     Xt is the first difference, and Xt.est is the estimate for Xt from
##'     EDM. Row numbers are time, and start at 1.
##'     Could maybe extend to have multiple Xt.est for the different
##'     embedding dimensions. TODO
##' @param pdf.filename filename (including .pdf) to save the movie as if saving
##'  as .pdf. Ignored if open.pdf is FALSE.
##' @param rhoForE values of rho corrsponding to each E in Evec
##' @param Evec vector of E values used to construct columns in Nx.lags called
##'     XtPredEeq1, XtPredEeq2 etc. (Xt predicted with E=1, etc.). If NULL
##'     then no predicted vs observed plot is drawn.
##' @param Ecols colour coding for E points in predicted vs observed plot
##' @param includeTimeSeries TRUE to include the time series at the top (FALSE
##'     hasn't  been tested yet TODO)
##' @param final_plot What the final plot should be. Currently only works for
##'   `compare` to compare Andy's predictions with rEDM's. Will then extend to
##'   compare different E values specified in `Evec` with `final_plot` being "E_vary".
##' @param start first time value (row) to use when plotting - not fully
##'     implemented as the colours aren't correct. Gives warning if != 1. TODO
##' @param end last time value (row) to use when plotting (so will have
##'     end-start points plotted) TODO (end-start-1??).
##' @param max_time maximum time value to plot (will be different from `end` if
##'     making a movie using Rmarkdown, since for that `end` will get iterated)
##' @param only.final.plot if TRUE then only plot the final plot
##' @param cobwebbing whether to join up points in 2d lag plots via cobwebbing;
##'    if FALSE then joins consecutive points (TODO?? - should be TRUE?).
##' @param late.col colour for the late points (the last 'late.num' points)
  #     for each frame
##' @param early.col colour for the early points
##' @param early.col.lines colour for the early lines
##' @param late.num how many points to have colour 'late.col'
##' @param pt.type usual 'type' option (but don't use "o")
##' @param x.lab x-label for axes for 3d lagged plots - for 2d
##'        lagged plots z.lab is the y label and y.lab is the x label.
##' @param y.lab y-label for axes for 3d lagged plots - for 2d
##'        lagged plots z.lab is the y label and y.lab is the x label.
##' @param z.lab y-label for axes for 3d lagged plots - for 2d
##'        lagged plots z.lab is the y label and y.lab is the x label.
##' @param figheight height of .pdf, but ignored if open.pdf == FALSE.
##' @param figwidth width of .pdf, but ignored if open.pdf == FALSE.
##' @param open.pdf whether or not to open a new .pdf file for the figure. Use
##'       FALSE if .pdf is being opened automatically in knitr (then figheight
##'       and figwidth get ignored, and should be set in <<..>>=.
##' @param ... Arguments to be passed to methods
##' @return Creates movie as a .pdf then can then be easily input into Sweave file
##'     using the animategraphics latex pacakge. TODO make .gif option also for vignette.
##' @export
##' @author Andrew Edwards
plotPanelMovie.df2 = function(Nx.lags = Nx_lags_orig,
                              pdf.filename = "Nx_lags_orig_movie.pdf",
                              rhoForE = NULL,
                              Evec = NULL,
                              Ecols = NA,
                              includeTimeSeries = TRUE,
                              final_plot = "compare",
                              start = 1,
                              end = nrow(Nx.lags),
                              max_time = NULL,
                              only.final.plot = FALSE,
                              cobwebbing = TRUE,
                              late.col = "red",
                              early.col = "black",
                              early.col.lines = "lightgrey",
                              late.num = 3,
                              pt.type = "p",
                              x.lab = expression("x"[t-2]),
                              y.lab = expression("x"[t-1]),
                              z.lab = expression("x"[t]),
                              figheight = 9.8,
                              figwidth = 6.84,
                              open.pdf = TRUE,
                              ...)
  {
#  max.time = nrow(Nx.lags.df)
#  time = 1:max.time
  if(start != 1) warning("Not properly implemented yet for start > 1;
                                colours won't work.")

  Nx.lags.use = Nx.lags  # not sure is needed

  # Axes ranges:
  if(is.null(max_time)){ max_time = end }
  t.axis.range = c(start, max_time)

  Nt.max.abs = max( abs( range(Nx.lags.use[start:max_time, "Xt"],
                                 na.rm=TRUE) ) )
  Nt.axes.range = c(0, Nt.max.abs*1.04)         # Expand else points can hit edge

  Xt.max.abs = max( abs( range(Nx.lags.use[start:max_time, "Xt"], na.rm=TRUE) ) )
  Xt.axes.range = c(-Xt.max.abs, Xt.max.abs)    # Make axes symmetric, though
                                                #  axs="i" doesn't work for 3d
  if(open.pdf) pdf(pdf.filename, height = figheight, width = figwidth)

  first.fig = ifelse( only.final.plot, end, start)
                                                # if only doing final plot then
                                                #  make just one figure in loop:
  for(iii in first.fig:end)    # Use iii here, which is now the row number
    {
      # Loops: iii=1 will only have pred vs obs,
      #        iii=2 will only have pred vs obs and 2d plots
      #        iii=3 onwards will have all plots
      ifelse(includeTimeSeries,
             par(mfrow = c(3,2)),
             par(mfrow = c(2,2)))
#      par.mai = c(0.7, 0.8, 0.1, 0.1)
#      par(mai = par.mai)# Usually affects all figures but scatterplot3d resets
                        #  so have to set this again later. Actually it affects
                        #  mar
      # Set mar, the numbers of lines of margins, default c(5, 4, 4, 2) + 0.1.
      par.mar.ts = c(3, 3, 1, 1)         # For time series
      par.mar.phase = c(3, 0, 1, 0)      # For phase plots (3d sets it anyway)

      par.mgp = c(1.5, 0.5, 0)
      par("mgp" = par.mgp) # first val sets axis title dist (mex units)
                           #  to axes, second sets labels


      i = iii-1       # to not have to change indexing from plotLag2d()

      # Colour vector for all plots - it will correspond
      #  to the last late.num times, not points (since different plots have
      #  different numbers of points). So these now correspond to times from
      #  start:iii , and so each needs to have length iii-start+1 (maybe not
      #  lines all get used):
      col.plot = c( rep(early.col, max(c(0, iii-late.num-start+1))),
                    rep(late.col, min(c(iii, late.num))) )   # colours of points
      col.plot.lines = col.plot                            # colours of lines
      col.plot.lines[col.plot.lines == early.col] = early.col.lines
      pch.plot = (col.plot == early.col) * 1 + (col.plot == late.col) * 16
                                         # filled circles for latest
      pch.plot[length(pch.plot)] = 8     # latest one a star

      if(includeTimeSeries)
        {
          par(pty = "m")                 # maximal plotting region, not square
                                         #  like for phase plots
          par(mar = par.mar.ts)
          # N_t vs t:
          plot(0, 0,
               xlab = expression("Time, t"),
               ylab = expression("N"[t]),
               xlim = t.axis.range,
               ylim = Nt.axes.range,
               type = "n",
               main = paste0("Time t=", iii))       # empty plot
          if(iii > 1.5)
            {
              segments(start:(iii-1),
                       dplyr::pull(Nx.lags.use[start:(iii-1), "Nt"]),
                       (start+1):iii,
                       dplyr::pull(Nx.lags.use[(start+1):iii, "Nt"]),
                       col = col.plot.lines) # lines() will not use vector col
            }
          points(start:iii,
                 dplyr::pull(Nx.lags.use[start:iii, "Nt"]),
                 type = pt.type,
                 pch = pch.plot,
                 col = col.plot)
           # x_t vs t:
          XtLoc = -0.05 * max(t.axis.range)  # location to plot Xt on a vertical line,
                                             #  needs correcting if start>1
          plot(0, 0,
               xlab = expression("Time, t"),
               ylab = expression("x"[t]),
               xlim = c(XtLoc, max(t.axis.range)),
               ylim = Xt.axes.range,
               type = "n")                           # empty plot
          abline(v = 0.5*XtLoc, col="black")
          if(iii > 1.5)
            {
              segments(start:(iii-1),
                       dplyr::pull(Nx.lags.use[start:(iii-1), "Xt"]),
                       (start+1):iii,
                       dplyr::pull(Nx.lags.use[(start+1):iii, "Xt"]),
                       col = col.plot.lines) # lines() will not use vector col
            }
          points(start:iii,
                 dplyr::pull(Nx.lags.use[start:iii, "Xt"]),
                 type = pt.type,
                 pch = pch.plot,
                 col = col.plot)
                                                     # '1d phase plot':
          points(rep(XtLoc, iii-start+1),
                 dplyr::pull(Nx.lags.use[start:iii, "Xt"]),
                 type = pt.type, pch = pch.plot,
                 col = col.plot)

        }                                            # end if(includeTimeSeries)

      par(pty="s")             # set following plot types to be square
                               #  (without this the axes don't seem to
                               #  be the same, even with the settings below)
      par(mar = par.mar.phase) # margins
      # N_t vs N{t-1}:
      # Empty plot to get started, that's it for iii=0:
      plot(0, 0,
           xlab = expression("N"[t-1]),
           ylab = expression("N"[t]),
           xlim = Nt.axes.range,
           ylim = Nt.axes.range,
           type = "n")
      if(cobwebbing) abline(0, 1, col="darkgrey")
      # Draw lines first so they get overdrawn by points
      if(iii > 2.5)
        {
          if(cobwebbing)
            {
              # Do lines for cobwebbing
              Nvals = rep(dplyr::pull(Nx.lags.use[start:iii, "Nt"]),
                          each = 2)
              Nvals = Nvals[-1]
              Nvals = Nvals[-length(Nvals)]
              len = length(Nvals)
              col.cobweb.lines = rep(early.col.lines, len)
              col.cobweb.lines[(max( (len - 2*late.num + 1), 1)):len] = late.col
              segments(Nvals[1:(len-2)],
                       Nvals[2:(len-1)],
                       Nvals[2:(len-1)],
                       Nvals[3:len],
                       col = col.cobweb.lines)
            } else
            {
              # Join each point to the next
              segments(dplyr::pull(Nx.lags.use[start:(iii-1), "Ntmin1"]),
                       dplyr::pull(Nx.lags.use[start:(iii-1), "Nt"]),
                       dplyr::pull(Nx.lags.use[(start+1):iii, "Ntmin1"]),
                       dplyr::pull(Nx.lags.use[(start+1):iii, "Nt"]),
                       col = col.plot.lines) # lines() will not use vector of col
            }
        }
      if(iii > 1.5)
        {
          points(dplyr::pull(Nx.lags.use[start:iii, "Ntmin1"]),
               dplyr::pull(Nx.lags.use[start:iii, "Nt"]),
               type = pt.type, pch = pch.plot,
               col = col.plot)          # start row has NA's, gets ignored
        }
      # legend("topright", legend=paste0("Time t=", iii), box.col = "white",
      #        inset = 0.01)  # inset to stop white overwriting outer box


      # x_t vs x{t-1}:

      # Empty plot to get started:
      plot(0, 0,
           xlab = y.lab, ylab = z.lab,
           xlim = Xt.axes.range,
           ylim = Xt.axes.range,
           type = "n")
      if(cobwebbing) abline(0, 1, col="darkgrey")
      if(iii > 2.5)
        {
          if(cobwebbing)
            {
               xvals = rep( dplyr::pull(Nx.lags.use[start:iii, "Xt"]), each = 2)
               xvals = xvals[-1]
               xvals = xvals[-length(xvals)]
               lenx = length(xvals)
               segments(xvals[1:(lenx-2)],
                        xvals[2:(lenx-1)],
                        xvals[2:(lenx-1)],
                        xvals[3:lenx],
                        col = col.cobweb.lines)
            } else
            {  # Just join consecutive points with lines
               segments(dplyr::pull(Nx.lags.use[start:(iii-1), "Xtmin1"]),
                        dplyr::pull(Nx.lags.use[start:(iii-1), "Xt"]),
                        dplyr::pull(Nx.lags.use[(start+1):iii, "Xtmin1"]),
                        dplyr::pull(Nx.lags.use[(start+1):iii, "Xt"]),
                        col = col.plot.lines)
            }
        }
      if(iii > 1.5)
        {
          points(dplyr::pull(Nx.lags.use[start:iii, "Xtmin1"]),
                 dplyr::pull(Nx.lags.use[start:iii, "Xt"]),
                 type = pt.type, pch = pch.plot,
                 col = col.plot)           # start row has NA's, get ignored
        }
      # legend("topleft", legend=paste("Time", iii), border = NULL)

      # 3d plot of x_t vs x_{t-1} vs x_{t-2}:
      # Empty plot to get started
      par.mgp.3d = c(3, 10, 0)
      par(mgp = par.mgp.3d)
      par(mai = c(0.1, 0.1, 0.1, 0.1)) # scat..3d resets mar, think mai still
      # has an effect
      par.mar.3d = c(3, 0, 0, 0)
      scat = scatterplot3d::scatterplot3d(0,
                                          0,
                                          0,
                                          xlab = x.lab,
                                          ylab = y.lab,
                                          zlab = z.lab,
                                          xlim = Xt.axes.range,
                                          ylim = Xt.axes.range,
                                          zlim = Xt.axes.range,
                                          type = "n",
                                          box = FALSE,
                                          angle = 40 + iii,
                                          mar = par.mar.3d)
                                              # mar=c(5,3,4,3)+0.1 is default,set
                                              #  within scatterplot3d
                                              # Add axes so can see origin:
      actual.axes.ranges = gets3dusr(scat)    # Get the actual values used
      scat$points3d(actual.axes.ranges[1:2], c(0,0), c(0,0), type = "l",
                    col = "lightgrey")
      scat$points3d(c(0,0), actual.axes.ranges[3:4], c(0,0), type = "l",
                    col = "lightgrey")
      scat$points3d(c(0,0), c(0,0), actual.axes.ranges[5:6], type = "l",
                    col = "lightgrey")
      # Obtain x-y co-ords of points for segments:
      proj.pts = scat$xyz.convert(
                                  dplyr::pull(Nx.lags.use[start:iii, "Xtmin2"]),
                                  dplyr::pull(Nx.lags.use[start:iii, "Xtmin1"]),
                                  dplyr::pull(Nx.lags.use[start:iii, "Xt"]) )
      if(iii > 3.5)
        {   # Think the indexing will now be 1:(iii-start), need start value also
            segments(proj.pts$x[1:(iii-start)], proj.pts$y[1:(iii-start)],
                    proj.pts$x[2:(iii-start+1)], proj.pts$y[2:(iii-start+1)],
                    col = col.plot.lines) # lines() will not use vector
                                          #  of col
        }
      # The points
      if(iii > 2.5)
        {
          scat$points3d(dplyr::pull(Nx.lags.use[start:iii, "Xtmin2"]),
                        dplyr::pull(Nx.lags.use[start:iii, "Xtmin1"]),
                        dplyr::pull(Nx.lags.use[start:iii, "Xt"]),
                        type = pt.type, pch = pch.plot,
                        col = col.plot)
        }


#      par(scat$par.mar)    # should do the same as:
      par(mar = par.mar.phase)   # scatterplot3d changes mar
      par(mgp = par.mgp)   # back to usual for 2d figures
      if(final_plot == "compare")
      {
        compare_cols = c("blue", "red")    # compare rEDM pred then Andy's.
        compare_pch = c(20, 1)
        pred_obs_max_abs = max( abs( range( c(Nx.lags$Xt,
                                              Nx.lags$rEDM.pred,
                                              Nx.lags$my.pred),
                                           na.rm = TRUE)
                                    )
                               )
        pred_obs_axes = c(-pred_obs_max_abs, pred_obs_max_abs)
        #  all.pred = dplyr::select(Nx.lags.use, starts_with("XtPredEeq"))
        #  pred.max.abs = max( abs( range(all.pred, na.rm=TRUE) ) )
        #  pred.max.abs = max(pred.max.abs, Xt.max.abs)  # Latter is observed
        #  predObs.axes.range = c(-pred.max.abs, pred.max.abs)

        plot(0,
             0,
             xlab = expression("Observation of x"[t]),
             ylab = expression("Prediction of x"[t]),
             xlim = pred_obs_axes,
             ylim = pred_obs_axes,
             asp = 1,
             type = "n")
        abline(0, 1, col="grey")
        leg = vector()
        points(dplyr::select(Nx.lags[start:iii,],
                             Xt,
                             rEDM.pred),
               pch = compare_pch[1],
               col = compare_cols[1])
        points(dplyr::select(Nx.lags[start:iii,],
                             Xt,
                             my.pred),
               pch = compare_pch[2],
               col = compare_cols[2])
        legend("topleft",
               pch = compare_pch,
               leg = c("rEDM pred", "Andy pred"),
               col = compare_cols)
      }

      # Predictions vs observations for E values in Evec
      if(final_plot == "E_vary")
        if(!is.null(Evec))
        {
          all.pred = dplyr::select(Nx.lags.use,
                                   dplyr::starts_with("XtPredEeq"))
          pred.max.abs = max( abs( range(all.pred, na.rm=TRUE) ) )
          pred.max.abs = max(pred.max.abs, Xt.max.abs)  # Latter is observed
          predObs.axes.range = c(-pred.max.abs, pred.max.abs)

          plot(0, 0,
               xlab = expression("Observation of x"[t]),
               ylab = expression("Prediction of x"[t]),
               xlim = predObs.axes.range,
               ylim = predObs.axes.range,
               asp = 1,
               type = "n")
          abline(0, 1, col="grey")
          leg = vector()
          for(j in 1:length(Evec))
            {
               points(dplyr::select(Nx.lags[start:iii,], Xt, paste0("XtPredEeq", j)),
                      pch=pch.plot, col=Ecols[j])
               leg = c(leg,
                       paste0("E=", j, ", rho=", round(rhoForE[j], 2)))
            }
          legend("topleft", pch=c(20, 20, 20),
                 leg, col=Ecols, cex=0.7)
        }                           # if(Evec != 0)
  }                                # for(iii in start:end)
  if(open.pdf) grDevices::dev.off() # close pdf device
}


##' Simplex prediction only for `E=2`, only used to run `data-raw/Nx_lags_orig.R`
##'
##' Simple early code to estimate \eqn{X(t^*+1)} and its variance, for
##'  all valid values of \eqn{t^*} for a `tbl_df` `Nx.lags`. Only for embedding dimension
##'  \eqn{E=2}. Was original older code, keeping in package to demonstrate what
##'  was done, and that this code is independent to the `pbsEDM()`
##'  function. Probably no need to make this code consistent with final
##'  notation, as this code is not mean to be used again.
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
