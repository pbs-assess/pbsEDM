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
##'   s3d <- scatterplot3d::scatterplot3d(rnorm(5), rnorm(5), rnorm(5))
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


##' Plot data and results of an object of class `pbsEDM`
##'
##' Only working if nt_observed values are in there (needs another switch if
##' not). And very likely only works for univariate time series for now, not
##' pred-prey for example.
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param portrait if TRUE then plots panels in portrait mode for manuscripts/Rmarkdown (3x2), false is
##'   landscape (2x3, filled columnwise) for presentations.
##' @param ... additiontal arguments to be passed onto other plotting functions,
##'   in particular `last.time.to.plot` to plot only up to that time step (to
##'   loop through in a movie), and late.num to plot the final number of time
##'   steps in a different colour
##' @return
##' @export
##' @author Andrew Edwards
plot.pbsEDM = function(obj,
                       portrait = TRUE,
                       ...){
  stopifnot(attr(obj, "class") == "pbsEDM")

  ifelse(portrait,
         par(mfrow = c(3,2)),
         par(mfcol = c(2,3)))

  if(!exists("last.time.to.plot")) last.time.to.plot <- length(obj$xt_observed)
                                            # though last will have NA # not
                                            # incorporated fully yet
  plot_observed(obj,
                ...)
#                end = last.time.to.plot)

  plot_phase_2d(values = obj$nt_observed,
                X.or.N = "N",
                ...)

  plot_phase_2d(values = obj$xt_observed,
                X.or.N = "X",
                ...)

  plot_phase_3d(obj,
                ...)
}

##' Plot output from pbsEDM_Evec as six panel plot including final forecast v
##'  observed for each E
##'
##'  <description>
##'
##' @return Six panel plot
##' @export
##' @author Andrew Edwards
##' @param E_res List of `pbsEDM` objects as output from `pbsEDM_Evec()`
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(Nx_lags_orig$Nt)
##'   plot_pbsEDM_Evec(aa)
##' }
plot_pbsEDM_Evec <- function(E_res,
                             ...){
  # The data are the same for each value of E, and so use one for the first five panels:
  plot.pbsEDM(E_res[[1]],
              ...)

  plot_pred_obs(E_res,
                ...)
}

##' Plot rho versus E for output from pbsEDM_Evec

##'  <description>
##'
##' @param E_res List of `pbsEDM` objects as output from `pbsEDM_Evec()`
##' @param ...
##' @return plot of rho versus E
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(Nx_lags_orig$Nt)
##'   plot_rho_Evec(aa)
##' }
plot_rho_Evec <- function(E_res,
                          ...){
  rho <- vector()
  E <- vector()
  for(i in 1:length(E_res)){
    rho[i] <- E_res[[i]]$results$rho
    E[i] <- E_res[[i]]$results$E
  }

  plot(E,
       rho,
       ylim = c(0,1),
       xlab = "Embedding Dimension (E)",
       ylab = "Forecast Skill (rho)",
       ...
       )
  # see EDMsimulate/report/sockeye-sim-edm.rnw for adding other plots
}


##' Plot predictions of X[t] versus the observations for list of `pbsEDM` objects
##'
##' <description>
##'
##' @return
##' @export
##' @author Andrew Edwards
##' @param E_res A list of `pbsEDM` objects
##' @param E_components How many of the first E_res components to show
##' @param E_cols Vector of colours, one for each value of E
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(Nx_lags_orig$Nt)
##'   plot_pred_obs(aa)
##' }
plot_pred_obs <- function(E_res,
                          E_components = 5,
                          E_cols = c("orange",
                                     "blue",
                                     "green",
                                     "red",
                                     "black"),
                          last.time.to.plot = NULL,
                          ...
                          ){

  stopifnot("Need more distinct colours in E_cols"=
              length(E_cols) <= E_components)

  if(is.null(last.time.to.plot)) last.time.to.plot <- length(E_res[[1]]$xt_observed)

  # Determine range for axes
  obs.max.abs = max( abs( range(E_res[[1]]$xt_observed,
                                na.rm=TRUE) ) ) # max abs observed value
  forecast.max.abs = 0         # max abs forecast value
  for(j in 1:E_components){
    forecast.max.abs = max(c(forecast.max.abs,
                             abs( range( range(E_res[[j]]$xt_forecast,
                                               na.rm=TRUE) ) ) ) )
  }

  max.abs = max(obs.max.abs, forecast.max.abs)
  axes.range = c(- max.abs, max.abs)

  plot(0, 0,
       xlab = expression("Observation of X"[t]),
       ylab = expression("Prediction of X"[t]),
       xlim = axes.range,
       ylim = axes.range,
       asp = 1,
       type = "n")
  abline(0, 1, col="grey")
  leg = vector()

  for(j in 1:E_components){
    points(E_res[[j]]$xt_observed[1:(last.time.to.plot-1)],
           E_res[[j]]$xt_forecast[1:(last.time.to.plot-1)],
           pch = 16, # could note the final one differently (but not a star,
                     # since that's last N[t] not X[t])
           col = E_cols[j])
    leg = c(leg,
            paste0("E=",
                   E_res[[j]]$results$E,
                   ", rho=",
                   round(E_res[[j]]$results$rho,
                         2)))
  }

  legend("topleft",
         pch=c(20, 20, 20),
         leg,
         col=E_cols,
         cex=0.7)
}


##' Plot values of X(t), X(t-1), X(t-2) in one of various ways
##'
##' Adapting from code within plotPanelMovie.df2. Assume raw data are one-dimensional
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param dim number of dimensions to plot: 1 for time series, 2 for X(t) vs X(t-2),
##'   3 for X(t-2) v X(t-1) v X(t) - not now TODO
##' @param late.col
##' @param early.col
##' @param early.col.lines
##' @param late.num
##' @param pt.type
##' @param x.lab
##' @param y.lab
##' @param z.lab
##' @param ...
##' @param last.time.to.plot last time value of N[t] to use when plotting, so
##'   final X[t] used is X[t-1] (since X[t] uses N[t+1])
##' @param max_time maximum time value for the time axis (will be different from `end` if
##'     making a movie using Rmarkdown, since for that `end` will get iterated)
##' @return
##' @export
##' @author Andrew Edwards
plot_observed = function(obj,
                         dim = 1,    # TODO remove and put into wrapper function
                                     # if needed
                         last.time.to.plot = NULL,
                         max_time = NULL,
                         late.col = "red",
                         early.col = "black",
                         early.col.lines = "lightgrey",
                         late.num = 3,
                         pt.type = "p",
                         x.lab = expression("X"[t-2]),
                         y.lab = expression("X"[t-1]),
                         z.lab = expression("X"[t]),
                         ...){
  stopifnot(attr(obj, "class") == "pbsEDM")
  stopifnot(dim %in% 1:3)

  # First row of everything is t=1
  if(is.null(max_time)) max_time <- length(obj$xt_observed)
  if(is.null(last.time.to.plot)) last.time.to.plot <- max_time

  start = 1    # start time of plots, only works for 1
  t.axis.range = c(start, max_time)

  if(!is.null(obj$nt_observed)){
    Nt.max.abs = max( abs( range(obj$nt_observed[start:max_time],
                                 na.rm=TRUE) ) )
    Nt.axes.range = c(0, Nt.max.abs*1.04)    # Expand else points can hit edge
  }

  Xt.max.abs = max(abs( range(obj$xt_observed[start:max_time], na.rm=TRUE) ),
                   abs( range(obj$xt_forecast[start:max_time], na.rm=TRUE)) )
  Xt.axes.range = c(-Xt.max.abs, Xt.max.abs)
  # Set mar, the numbers of lines of margins, default c(5, 4, 4, 2) + 0.1.
  #  par.mar.ts = c(3, 3, 1, 1)         # For time series
  par.mar.phase = c(3, 0, 1, 0)      # For phase plots (3d sets it anyway)
  par.mar.3d = c(3, 0, 0, 0)

  par.mgp.3d = c(3, 10, 0)
  par.mgp = c(1.5, 0.5, 0)
  par("mgp" = par.mgp) # first val sets axis title dist (mex units)
                       #  to axes, second sets labels

  # iii = last.time.to.plot - 1    # depends on the type of plot, don't have value of
                   # Xt[last.time.to.plot]. This is how many points to plot.

  # Colour vector for all plots - it will correspond
  #  to the last late.num times, not points (since different plots have
  #  different numbers of points).
  # Redoing 16/6/20 - may be different to my old way. So plot the data up to
  # time last.time.to.plot, but in the sense of using N[last.time.to.plot] but
  # not X[last.time.to.plot] because that uses N[last.time.to.plot + 1]
  col.plot = c(rep(early.col,
                   max(c(0, last.time.to.plot - late.num))),
               rep(late.col,
                   min(c(last.time.to.plot, late.num))) )   # colours of points
  col.plot.lines = col.plot                            # colours of lines
  col.plot.lines[col.plot.lines == early.col] = early.col.lines
  pch.plot = (col.plot == early.col) * 1 + (col.plot == late.col) * 16
                                        # filled circles for latest
  pch.plot[length(pch.plot)] = 8     # latest one a star

  # Nt v t (if available) and Xt v t, with Xt points also shown on 1-d line
  if(dim == 1){
    if(!is.null(obj$nt_observed)){
    #  par(mfrow=c(1,2))       # Will have to generalise this

      plot_time_series(values = obj$nt_observed[start:last.time.to.plot],
                       X.or.N = "N",
                       par.mar.ts = c(3, 3, 1, 1),
                       t.axis.range = t.axis.range,
                       y.range = Nt.axes.range,
                       col.plot = col.plot,
                       col.plot.lines = col.plot.lines,
                       pch.plot = pch.plot,
                       last.time.to.plot = last.time.to.plot   # For title
                       )
    }

    plot_time_series(values = obj$xt_observed[start:(last.time.to.plot-1)], # don't want to use N[last.time.to.plot]
                     X.or.N = "X",
                     par.mar.ts = c(3, 3, 1, 1),
                     t.axis.range = t.axis.range,
                     y.range = Xt.axes.range,
                     col.plot = col.plot,
                     col.plot.lines = col.plot.lines,
                     pch.plot = pch.plot
                     )
  }
}

##' Plot the observed time series as either Nt or Xt
##'
##'  <description>
##'
##' @param values
##' @param X.or.N "N" if raw non-differenced data, "X" for differenced data
##' @param par.mar.ts
##' @param t.axis.range
##' @param iii time value to plot up to
##' @param y.range
##' @param col.plot
##' @param col.plot.lines
##' @param pch.plot
##' @param start
##' @param pt.type
##' @return
##' @export
##' @author Andrew Edwards
plot_time_series <- function(values,
                             X.or.N,
                             par.mar.ts,
                             t.axis.range = NULL,
                             last.time.to.plot = NA,
                             y.range,
                             col.plot,
                             col.plot.lines,
                             pch.plot,
                             start = 1,
                             pt.type = "p"
                             ){
  if(is.null(t.axis.range)) {
    t.axis.range <- c(1, length(values))
  }
  if(is.na(last.time.to.plot)) {
    last.time.to.plot = max(t.axis.range)
  }

  par(pty = "m")                 # maximal plotting region, not square
                                 #  like for phase plots
  par(mar = par.mar.ts)

  if(X.or.N == "N"){
    plot(0, 0,
         xlab = expression("Time, t"),
         ylab = expression("N"[t]),
         xlim = c(0, max(t.axis.range)),
         ylim = y.range,
         type = "n",                           # empty plot
         main = paste0("Time t=", last.time.to.plot))
  } else {
    XtLoc = -0.05 * max(t.axis.range)  # location to plot Xt on a vertical line,
    plot(0, 0,
         xlab = expression("Time, t"),
         ylab = expression("X"[t]),
         xlim = c(XtLoc, max(t.axis.range)),
         ylim = y.range,
         type = "n")                           # empty plot
    abline(v = 0.5*XtLoc, col="black")
  }

  iii = last.time.to.plot             # needs replacing when I have time
  if(iii > 1.5){
    segments(start:(iii-1),
             values[start:(iii-1)],    # dplyr::pull(Nx.lags.use[start:(iii-1), "Xt"]),
             (start+1):iii,
             values[(start+1):iii],
             col = col.plot.lines)      # lines() will not use vector col
  }

  points(start:iii,
         values[start:iii],
         type = pt.type,
         pch = pch.plot,
         col = col.plot)

  # '1d phase plot':
  if(X.or.N == "X"){
    points(rep(XtLoc, iii-start+1),
           values[start:iii],
           type = pt.type,
           pch = pch.plot,
           col = col.plot)
  }
}
##' @
##'  <description>
##'
##' @param values
##' @param X.or.N
##' @param par.mar.phase
##' @param axis.range
##' @param iii
##' @param pt.type
##' @param cobwebbing
##' @param late.col
##' @param early.col
##' @param early.col.lines
##' @param late.num
##' @param y.lab
##' @param z.lab
##' @param ... additiontal arguments,
##'   in particular `last.time.to.plot` to plot only up to that time step (to
##'   loop through in a movie).
##' @return
##' @export
##' @author Andrew Edwards
plot_phase_2d <- function(values,
                          X.or.N = "X",
                          par.mar.phase = c(3, 0, 1, 0),
                          axis.range = NA,
                          start = 1,
                          last.time.to.plot = NULL,
                          pt.type = "p",
                          cobwebbing = TRUE,
                          late.col = "red",
                          early.col = "black",
                          early.col.lines = "lightgrey",
                          late.num = 3,
                          y.lab = expression("x"[t-1]),
                          z.lab = expression("x"[t]),
                          ...
                          ){

  if(is.na(axis.range)) {
    axis.range <- c(min(0, min(values, na.rm = TRUE)),
                    max(values, na.rm = TRUE))   # for Nt or Xt
  }
  if(is.null(last.time.to.plot)) {
    last.time.to.plot = length(values)
  }

  # Now copying from plot_observed:
  col.plot = c(rep(early.col,
                   max(c(0, last.time.to.plot - late.num))),
               rep(late.col,
                   min(c(last.time.to.plot, late.num))) )   # colours of points
  col.plot.lines = col.plot                            # colours of lines
  col.plot.lines[col.plot.lines == early.col] = early.col.lines
  pch.plot = (col.plot == early.col) * 1 + (col.plot == late.col) * 16
                                        # filled circles for latest
  pch.plot[length(pch.plot)] = 8     # latest one a star

  # Copying from plotPanelMovie.df2

  par(pty="s")             # set following plot types to be square
                           #  (without this the axes don't seem to
                           #  be the same, even with the settings below)
  par(mar = par.mar.phase) # margins

  # N_t vs N{t-1}:
  # Empty plot to get started
  if(X.or.N == "N"){
    plot(0, 0,
         xlab = expression("N"[t-1]),
         ylab = expression("N"[t]),
         xlim = axis.range,
         ylim = axis.range,
         type = "n")
    values.to.plot <- values[start:last.time.to.plot]
  } else {
    plot(0, 0,
         xlab = y.lab,
         ylab = z.lab,
         xlim = axis.range,
         ylim = axis.range,
         type = "n")
    values.to.plot <- values[start:(last.time.to.plot-1)] # Not use X[last..]
  }

  if(cobwebbing) abline(0, 1, col="darkgrey")

  # Draw lines first so they get overdrawn by points
  if(last.time.to.plot > 2.5){
    if(cobwebbing){
      # Do lines for cobwebbing
      vals = rep(values.to.plot,
                 each = 2)
      vals = vals[-1]
      vals = vals[-length(vals)]
      len = length(vals)
      col.cobweb.lines = rep(early.col.lines, len)
      col.cobweb.lines[(max(len - 2*late.num + 1 + 2*(X.or.N == "X"),
                            1)):len] = late.col
      segments(vals[1:(len-2)],
               vals[2:(len-1)],
               vals[2:(len-1)],
               vals[3:len],
               col = col.cobweb.lines)
    } else {
      # Join each point to the next   (N(5), N(6)) to (N(6), N(7))
      #  Not checked.
      segments(
        # dplyr::pull(Nx.lags.use[start:(iii-1), "Ntmin1"]),
        # dplyr::pull(Nx.lags.use[start:(iii-1), "Nt"]),
        # dplyr::pull(Nx.lags.use[(start+1):iii, "Ntmin1"]),
        # dplyr::pull(Nx.lags.use[(start+1):iii, "Nt"]),
        # not fully tested:
        pbsLAG(values.to.plot)[-length(values.to.plot)],   # N(t-1)
        values.to.plot[-length(values.to.plot)],
        pbsLAG(values.to.plot)[-1],
        values.to.plot[-1],
        col = col.plot.lines) # lines() will not use vector of col
    }
  }
  if(last.time.to.plot > 1.5){
    points(pbsLAG(values.to.plot),
           values.to.plot,
           type = pt.type,
           pch = pch.plot,
           col = col.plot)          # start row has NA's, gets ignored
  }
      # legend("topright", legend=paste0("Time t=", iii), box.col = "white",
      #        inset = 0.01)  # inset to stop white overwriting outer box


      # x_t vs x{t-1}:

#      if(iii > 2.5)
#        {
#          if(cobwebbing)
#            {
#               xvals = rep( dplyr::pull(Nx.lags.use[start:iii, "Xt"]), each = 2)
#               xvals = xvals[-1]
#               xvals = xvals[-length(xvals)]
#               lenx = length(xvals)
#              segments(xvals[1:(lenx-2)],
#                        xvals[2:(lenx-1)],
#                       xvals[2:(lenx-1)],
#                        xvals[3:lenx],
#                        col = col.cobweb.lines)
#           } else
#           {  # Just join consecutive points with lines
#               segments(dplyr::pull(Nx.lags.use[start:(iii-1), "Xtmin1"]),
#                        dplyr::pull(Nx.lags.use[start:(iii-1), "Xt"]),
#                        dplyr::pull(Nx.lags.use[(start+1):iii, "Xtmin1"]),
#                        dplyr::pull(Nx.lags.use[(start+1):iii, "Xt"]),
#                        col = col.plot.lines)
#            }
#        }
#      if(iii > 1.5)
#        {
#          points(dplyr::pull(Nx.lags.use[start:iii, "Xtmin1"]),
#                 dplyr::pull(Nx.lags.use[start:iii, "Xt"]),
#                 type = pt.type, pch = pch.plot,
#                 col = col.plot)           # start row has NA's, get ignored
#        }
      # legend("topleft", legend=paste("Time", iii), border = NULL)

}

##'  3d plot of x_t vs x_{t-1} vs x_{t-2}:
##'
##' <description>
##'
##' @param obj `pbsEDM` object
##' @param par.mgp.3d
##' @param par.mai.3d
##' @param par.mar.3d
##' @param x.lab
##' @param y.lab
##' @param z.lab
##' @param axis.range
##' @param late.num
##' @param pt.type
##' @param late.col
##' @param early.col
##' @param early.col.lines
##' @param axes.col
##' @param par.mar.phase
##' @return
##' @export
##' @author Andrew Edwards
plot_phase_3d <- function(obj,
                          par.mgp.3d = c(3, 10, 0),
                          par.mai.3d = c(0.1, 0.1, 0.1, 0.1),
                          par.mar.3d = c(3, 0, 0, 0),
                          x.lab = expression("x"[t-2]),
                          y.lab = expression("x"[t-1]),
                          z.lab = expression("x"[t]),
                          last.time.to.plot = NULL,
                          axis.range = NA,
                          late.num = 3,
                          pt.type = "p",
                          late.col = "red",
                          early.col = "black",
                          early.col.lines = "lightgrey",
                          axes.col = "darkblue",
                          par.mar.phase = c(3, 0, 1, 0),  # to reset for normal figs
                          par.mgp = c(1.5, 0.5, 0),
                          ...
                          ){
  if(is.na(axis.range)) {
    axis.range <- c(min(0, min(obj$xt_observed, na.rm = TRUE)),
                    max(obj$xt_observed, na.rm = TRUE))
  }

  if(is.null(last.time.to.plot)) last.time.to.plot <- length(obj$xt_observed)

  start = 1     # currently only works for 1

  if(last.time.to.plot > 2){
    values.to.plot <- obj$xt_observed[start:(last.time.to.plot-1)]   # can't use N[last...])
  } else {
    values.to.plot <- NA          # nothing to plot
  }

  # Copied from plot_observed:
  col.plot = c(rep(early.col,
                   max(c(0, last.time.to.plot - late.num))),
               rep(late.col,
                   min(c(last.time.to.plot, late.num))) )   # colours of points
  col.plot.lines = col.plot                            # colours of lines
  col.plot.lines[col.plot.lines == early.col] = early.col.lines
  pch.plot = (col.plot == early.col) * 1 + (col.plot == late.col) * 16
                                        # filled circles for latest
  pch.plot[length(pch.plot)] = 8     # latest one a star

  # Empty plot to get started
  par(mgp = par.mgp.3d)
  par(mai = par.mai.3d)  # scat..3d resets mar, think mai still has an effect
  par.mar.3d = c(3, 0, 0, 0)
  scat = scatterplot3d::scatterplot3d(0,
                                      0,
                                      0,
                                      xlab = x.lab,
                                      ylab = y.lab,
                                      zlab = z.lab,
                                      xlim = axis.range,
                                      ylim = axis.range,
                                      zlim = axis.range,
                                      type = "n",
                                      box = FALSE,
                                      angle = 40 + last.time.to.plot,
                                      mar = par.mar.3d)
                                              # mar=c(5,3,4,3)+0.1 is default,set
                                              #  within scatterplot3d
                                              # Add axes so can see origin:
      actual.axes.ranges = gets3dusr(scat)    # Get the actual values used
      scat$points3d(actual.axes.ranges[1:2], c(0,0), c(0,0), type = "l",
                    col = axes.col)
      scat$points3d(c(0,0), actual.axes.ranges[3:4], c(0,0), type = "l",
                    col = axes.col)
      scat$points3d(c(0,0), c(0,0), actual.axes.ranges[5:6], type = "l",
                    col = axes.col)
      # Obtain x-y co-ords of points for segments:
#**      proj.pts = scat$xyz.convert(dplyr::pull(Nx.lags.use[start:iii,
#                                                            "Xtmin2"]),
#                                    dplyr::pull(Nx.lags.use[start:iii,
#                                                            "Xtmin1"]),
  #                                    dplyr::pull(Nx.lags.use[start:iii, "Xt"]) )
  # maybe this is okay with NA's:
#  proj.pts = scat$xyz.convert(pbsLAG(obj$xt_observed,
#                                     2)[start:iii],  # "Xtmin2", index wrong?
#                              pbsLAG(obj$xt_observed,
#                                     1)[start:iii],  # "Xtmin1"
  #                              pbsLAG(obj$xt_observed)[start:iii])  # "Xt"
  if(all(!is.na(values.to.plot))){
    proj.pts = scat$xyz.convert(pbsLAG(values.to.plot,
                                       2),  # "Xtmin2"
                                pbsLAG(values.to.plot,
                                       1),  # "Xtmin1"
                                values.to.plot)  # "Xt"
    if(last.time.to.plot > 3.5){
      segments(proj.pts$x[1:(last.time.to.plot - start)],
               proj.pts$y[1:(last.time.to.plot - start)],
               proj.pts$x[2:(last.time.to.plot - start + 1)],
               proj.pts$y[2:(last.time.to.plot - start + 1)],
               col = col.plot.lines) # lines() will not use vector
    }

    # The points
    if(last.time.to.plot > 2.5){
      scat$points3d(pbsLAG(values.to.plot,
                           2),  # "Xtmin2"
                    pbsLAG(values.to.plot,
                           1),  # "Xtmin1"
                    values.to.plot,  # "Xt"
                    type = pt.type,
                    pch = pch.plot,
                    col = col.plot)
    }
  }

  par(mar = par.mar.phase)   # scatterplot3d changes mar
  par(mgp = par.mgp)        # back to usual for 2d figures
}

##' 2-d phase plot of x(t) v x(t-1) with coloured points to explain EDM
##'
##' Highlights a point to be projected, its nearest `E+1` neighbours, and then
##' draw arrows to show where they go and so where the projection goes. Very
##' useful for understanding and checking what EDM is doing. Call this multiple
##' times using TODO to make a movie. Type `plot_explain_edm_movie` for example calls.
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param tstar the time index (x(t)) of the target point for which to make a
##'   projection from
##' @return
##' @export
##' @author Andrew Edwards
plot_explain_edm <- function(obj,
                             tstar,
                             x.lab = expression("X(t-1)"),
                             y.lab = expression("X(t)"),
                             main = "All the points in lagged space",
                             tstar.col = "blue",
                             tstar.pch = 1,
                             tstar.cex = 1,
                             psivec = NULL,
                             neigh.plot = FALSE,
                             neigh.proj = FALSE,
                             pred.plot = FALSE,
                             pred.rEDM = FALSE,
                             true.val = FALSE,
                             legend.plot = TRUE){
  par(pty="s")

  Xt <- obj$xt_lags[, "Nt_0"]      # Should change Nt in pbsEDM() as confusing
  Xtmin1 <- obj$xt_lags[, "Nt_1"]

  Xt.max.abs = max(abs( range( c(Xt,
                                 Xtmin1,
                                 obj$xt_forecast),
                              na.rm=TRUE)))
  Xt.axes.range = c(-Xt.max.abs, Xt.max.abs)

  # Plot all the points
  plot(Xtmin1,
       Xt,
       xlab = x.lab,
       ylab = y.lab,
       xlim = Xt.axes.range,
       ylim = Xt.axes.range,
       main = main)

  # Highlight x[tstar]
  points(Xtmin1[tstar],
         Xt[tstar],
         col = tstar.col,
         pch = tstar.pch,
         cex = tstar.cex)

  psi_vec <- obj$xt_nbr_index[tstar,]     #  vector of indices of points closest
                                        #  to x[tstar]

  # Highlight neighbours
  if(neigh.plot){
    points(Xtmin1[psi_vec],
           Xt[psi_vec],
           pch = 19,
           col = "red")
  }

  # Highlight projections of neighbours
  if(neigh.proj){
    points(Xtmin1[psi_vec + 1],
           Xt[psi_vec + 1],
           pch = 19,
           col = "purple",
           cex = 0.5)

    for(i in 1:length(psi_vec)){
      igraph:::igraph.Arrows(Xtmin1[psi_vec[i]],
                              Xt[psi_vec[i]],
                              Xtmin1[psi_vec[i] + 1],
                              Xt[psi_vec[i] + 1],
                              curve = 1,
                              size = 0.7,
                              h.lwd = 2,
                              sh.lwd = 2,
                              width = 1,
                              sh.col = "purple")
               # Example:
               # iArrows(0, 0, 4, 4, curve = 1, size = 0.7, h.lwd = 2,
      #  sh.lwd=2, width = 1, sh.col="red")
    }
  }


  if(pred.plot){
    points(Xt[tstar],
           obj$xt_forecast[tstar + 1],    # CHECK not tstar+1, think it's correct
           col = tstar.col,
           pch = 8)

    igraph:::igraph.Arrows(Xtmin1[tstar],
                           Xt[tstar],
                           Xt[tstar],
                           obj$xt_forecast[tstar + 1],
                           curve = 1,
                           size = 0.7,
                           h.lwd = 2,
                           sh.lwd = 2,
                           width = 1,
                           sh.col = "darkgrey")

    abline(h = Xt[psi_vec + 1],
           col = "lightgrey")
  }

  if(pred.rEDM){
    points(Xt[tstar],
           dplyr::pull(Nx_lags_orig[tstar+1, "rEDM.pred"]),
           col = "red",
           cex = 1.5,
           lwd = 2)
  }

  if(true.val){
    points(Xtmin1[tstar + 1],
           Xt[tstar + 1],
           cex = 1.5,
           lwd = 2,
           col = "darkgreen")
  }

  if(legend.plot){
    legend("bottomleft",
           pch=c(tstar.pch, 19, 8, 1, 1),
           leg=c("x(t*)",
                 "neighbours",
                 "x(t*+1) pred (wt avge)",
                 "rEDM pred",
                 "true x(t*+1)"),
           col=c(tstar.col,
                 "red",
                 tstar.col,
                 "red",
                 "darkgreen"),
           cex=0.85)

    legend("topleft",
           leg = paste0("t*=", tstar),
           col = "white",
           bty = "n")
  }
}

##' Create movie to explain EDM
##'
##' Use with gifski in .Rmd - see vignette.
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param ... needs to include tstar
##' @return
##' @export
##' @author Andrew Edwards
##' \donttest{
##'   aa <- pbsEDM_Evec(Nx_lags_orig$Nt)
##'   plot_explain_edm_movie(aa, tstar = 15)
##' }
plot_explain_edm_movie <- function(obj,
                                   ...){
#browser()
#  stopifnot("Need to specify tstar" =
#              exists("tstar"))   # think this fails as it can't be used in THIS function

  plot_explain_edm(obj,
                   tstar.col = "black",
                   ...)    # tstar should get carried through

  #plot_explain_edm(obj,
  #                 tstar.col = "blue",,
  #                 main = paste0(
  #                   "Remove vector x(tstar) from library"),
  #                   # where tstar=", tstar,
  #                 tstar.pch = 4,
  #                 tstar.cex = 2.5,
  #                 ...)

  # DON'T we want to remove tstar+1 also - whole point of my error finding!
# TODO ADD THAT IN HERE

  #plot_explain_edm(obj,
  #                 tstar.col = "white",   # to hide it
  #                 main = paste0(
  #                   "Want to predict where it goes, i.e. predict vec x(tstar+1)"),
  #                 tstar.cex = 2.5,
  #                 ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "Want to predict where vector x(tstar) will go ..."),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "... based on three (E+1) nearest neighbours"),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   neigh.plot = TRUE,
                   ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "See where neighbours go"),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   neigh.plot = TRUE,
                   neigh.proj = TRUE,
                   ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "and take weighted average of X(t) to be predicted value"),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   neigh.plot = TRUE,
                   neigh.proj = TRUE,
                   pred.plot = TRUE,
                   ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "which should agree with prediction from rEDM"),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   neigh.plot = TRUE,
                   neigh.proj = TRUE,
                   pred.plot = TRUE,
                   pred.rEDM = TRUE,
                   ...)

  plot_explain_edm(obj,
                   main = paste0(
                     "and is hopefully close to the true value"),
                   tstar.pch = 19,
                   tstar.cex = 1.2,
                   neigh.plot = TRUE,
                   neigh.proj = TRUE,
                   pred.plot = TRUE,
                   pred.rEDM = TRUE,
                   true.val = TRUE,
                   ...)
}
