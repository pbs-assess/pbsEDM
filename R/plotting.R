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
##' @param ...
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
      par.mar.3d = c(3, 0, 0, 0)

      par.mgp.3d = c(3, 10, 0)
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
      par(mgp = par.mgp.3d)
      par(mai = c(0.1, 0.1, 0.1, 0.1)) # scat..3d resets mar, think mai still has an effect
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

#      for(j in 1:nrows(Nx.lags))
#            {
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

#               leg = c(leg,
#                       paste0("rEDM pred")E=", j, ", rho=", round(rhoForE[j], 2)
#            }
      legend("topleft",
             pch = compare_pch,
             leg = c("rEDM pred", "Andy pred"),
             col = compare_cols)
#             cex = 0.7)
      }
      # Predictions vs observations for E values in Evec
      if(final_plot == "E_vary")
        if(!is.null(Evec))
        {
          all.pred = dplyr::select(Nx.lags.use, starts_with("XtPredEeq"))
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
  if(open.pdf) dev.off() # close pdf device
}


# Obtain axes ranges, from Uwe Ligges, adapted from
#
##' Obtain axes ranges (adapted from function by Uwe Ligges)
##'
##' Obtain axes ranges for scatterplot3d plot (and maybe does others), adapted
##' from http://r.789695.n4.nabble.com/Axis-Limits-in-Scatterplot3d-td892026.html<description>
##'
##' @param s3dobject scatterplot3d plot presumably
##' @return Ranges of x, y and z axes
##' @export
##' @author Andrew Edwards (adapted from Uwe Ligges' function)
##' @examples
##'   s3d <- scatterplot3d(rnorm(5), rnorm(5), rnorm(5))
##'   gets3dusr(s3d)
gets3dusr = function(s3dobject)
  {
  es3d = environment(s3dobject$points3d)
  g3d  = function(name) get(name, envir=es3d)
  c(c(g3d("x.min"),
      g3d("x.max")) * g3d("x.scal"),
    c(0,
      g3d("y.max")) * g3d("y.scal") + g3d("y.add"),
    c(g3d("z.min"),
      g3d("z.max")) * g3d("z.scal"))
}
