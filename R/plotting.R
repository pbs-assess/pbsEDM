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
##' \donttest{
##'   s3d <- scatterplot3d::scatterplot3d(rnorm(5), rnorm(5), rnorm(5))
##'   gets3dusr(s3d)
##' }
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


##' Plot data and results of an object of class `pbsEDM` as time series and
##' phase plots
##'
##' Plots time series of `N_t` and `Y_t`, 2d phase plots of `N_t` vs `N_{t-1}`
##' and `Y_t` vs `Y_{t-1}`, and 3d phase plot of `Y_t` vs `Y_{t-1}` vs `Y_{t-2}`.
##' See vignette "Analyse a simple time series".
##' Only working if `N$N_t` values are in there (needs another switch if
##' not). And very likely only works for univariate time series for now, not
##' pred-prey for example. TODO, update if necessary. Works on `NY_lags_example`
##' but likely need to generalise.
##'
##' @param x list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param portrait if TRUE then plots panels in portrait mode for manuscripts/Rmarkdown (3x2), false is
##'   landscape (2x3, filled columnwise) for presentations.
##' @param ... additional arguments to be passed onto other plotting functions,
##'   in particular `last.time.to.plot` to plot only up to that time step (to
##'   loop through in a movie), and late.num to plot the final number of time
##'   steps in a different colour
##' @return Five panels in 2x3 or 3x2 format, in current plot environment
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'  aa <- pbsEDM(NY_lags_example,
##'               lags = list(N_t = 0:1),
##'               first_difference = TRUE)
##'  plot(aa)
##' }
plot.pbsEDM = function(x,
                       portrait = TRUE,
                       ...){
  stopifnot(attr(x, "class") == "pbsEDM")

  ifelse(portrait,
         par(mfrow = c(3,2)),
         par(mfcol = c(2,3)))

  # Commenting as isn't passed on anyway
  #if(!exists("last.time.to.plot")) last.time.to.plot <- length(x$X_observed)
  #                                          # though last will have NA # not
  #                                          # incorporated fully yet

  plot_time_series(values = x$N,  # $N_t,
                                        # TODO: why might as.vector() be needed , xt is okay
                   X.or.N = "N",
                   ...)

  plot_time_series(values = x$X_observed,
                   X.or.N = "X",
                   ...)

  plot_phase_2d(values = x$N,    # $N_t,
                X.or.N = "N",
                ...)

  plot_phase_2d(values = x$X_observed,
                X.or.N = "X",
                ...)

  plot_phase_3d(x,
                ...)
  invisible()
}

##' Plot output from `pbsEDM_Evec()` as six panel plot, with data showns the
##'  same as for `plot.pbsEDM()` and sixth panel with predictions vs observation
##'  for each `E`
##'
##' Plots time series of `N_t` and `Y_t`, 2d phase plots of `N_t` vs `N_{t-1}`
##' and `Y_t` vs `Y_{t-1}`, and 3d phase plot of `Y_t` vs `Y_{t-1}` vs
##' `Y_{t-2}`. Sixth panel shows predictions vs observations for different values
##' of `E`. See vignette "Analyse a simple time series" for full details.
##'
##' @param E_res List of `pbsEDM` objects as output from `pbsEDM_Evec()`
##' @param ... extra arguments to `plot.pbsEDM()`
##' @return Six-panel figure in current plot environment
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa_Evec <- pbsEDM_Evec(NY_lags_example$N_t)
##'   plot_pbsEDM_Evec(aa_Evec)
##'
##' # For manuscript figure, after running `analyse_simple_time_series.Rmd` vignette run:
##' postscript("six_panels.eps",
##'             height = 5.36,
##'             width = 9,
##'             horizontal=FALSE,
##'             paper="special")
##' for(iiii in 1:length(NY_lags_example$N_t)){
##'     plot_pbsEDM_Evec(E_results,
##'     last.time.to.plot = iiii,
##'     portrait = FALSE)
##'   }
##' dev.off()
##' }
plot_pbsEDM_Evec <- function(E_res,
                             ...){
  # The data are the same for each value of E, and so use one for the first five panels:
  plot.pbsEDM(E_res[[1]],
              ...)

  plot_pred_obs(E_res)
  invisible()
}

##' Plot rho versus `E` for output from `pbsEDM_Evec()`
##'
##' Plot correlation coefficient against `E`.
##'
##' @param E_res List of `pbsEDM` objects as output from `pbsEDM_Evec()`
##' @param ... Extra arguments to pass to `plot()`
##' @return plot of rho versus E
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(NY_lags_example$N_t)
##'   plot_rho_Evec(aa)
##' }
plot_rho_Evec <- function(E_res,
                          ...){
  rho <- vector()
  E <- vector()
  for(i in 1:length(E_res)){
    rho[i] <- E_res[[i]]$results$X_rho
    E[i] <- E_res[[i]]$results$E
  }

  plot(E,
       rho,
       ylim = c(0,1),
       xlab = "Embedding Dimension (E)",
       ylab = "Forecast Skill (rho)",
       ...
       )
  # TODO see EDMsimulate/report/sockeye-sim-edm.rnw for adding other plots
  invisible()
}


##' Plot predictions versus observed values of `Y_t` for list of `pbsEDM`
##' objects, for various values of `E`
##'
##' Plot the predicted versus observed values of `Y_t` for various values
##' of `E`, using the output (a list of `pbsEDM` objects) from `pbsEDM_Evec()`.
##'
##' @param E_res A list of `pbsEDM` objects
##' @param E_components How many of the first E_res components to show
##' @param E_cols Vector of colours, one for each value of E
##' @param last.time.to.plot Last time value of `N_t` to use when plotting.
##' @return Figure in the current plot enivironment
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM_Evec(NY_lags_example$N_t)
##'   plot_pred_obs(aa)
##' }
plot_pred_obs <- function(E_res,
                          E_components = 5,
                          E_cols = c("orange",
                                     "blue",
                                     "green",
                                     "red",
                                     "black"),
                          last.time.to.plot = NULL
                          ){

  stopifnot("Need more distinct colours in E_cols"=
              length(E_cols) <= E_components)

  if(is.null(last.time.to.plot)) last.time.to.plot <- length(E_res[[1]]$X_observed)

  # Determine range for axes
  obs.max.abs = max( abs( range(E_res[[1]]$X_observed,
                                na.rm=TRUE) ) ) # max abs observed value
  forecast.max.abs = 0         # max abs forecast value
  for(j in 1:E_components){
    forecast.max.abs = max(c(forecast.max.abs,
                             abs( range( range(E_res[[j]]$X_forecast,
                                               na.rm=TRUE) ) ) ) )
  }

  max.abs = max(obs.max.abs, forecast.max.abs)
  axes.range = c(- max.abs, max.abs)

  plot(0, 0,
       xlab = expression("Observation of Y"[t]),
       ylab = expression("Prediction of Y"[t]),
       xlim = axes.range,
       ylim = axes.range,
       asp = 1,
       type = "n")
  abline(0, 1, col="grey")
  leg = vector()

  for(j in 1:E_components){
    points(E_res[[j]]$X_observed[1:(last.time.to.plot-1)],
           E_res[[j]]$X_forecast[1:(last.time.to.plot-1)],
           pch = 16, # could note the final one differently (but not a star,
                     # since that's last N[t] not Y[t])
           col = E_cols[j])
    leg = c(leg,
            paste0("E=",
                   E_res[[j]]$results$E,
                   ", rho=",
                   round(E_res[[j]]$results$X_rho,
                         2)))
  }

  legend("topleft",
         pch=c(20, 20, 20),
         leg,
         col=E_cols,
         cex=0.7)
  invisible()
}

##' Plot the observed time series as either `N_t` or `Y_t`
##'
##' First value must be `t=1`. For non-differenced values `N_t`, shows the time
##' series with the final values in a different colour, and a title showing the
##' final time step. For first-differenced values `Y_t`, shows the time series
##' of those, plus a one-dimensional phase plot.
##'
##' @param values vector of values to be plotted
##' @param X.or.N "N" if raw non-differenced data, "X" for differenced data
##' @param par.mar.ts `par(mar)` values
##' @param max_time maximum time value for the time axis
##' @param t.axis.range range of time axis
##' @param last.time.to.plot last time value of N[t] to use when plotting, so
##'   final Y[t] used will be Y[t-1] (since Y[t] uses N[t+1])
##' @param late.num final number of `N[t]` time steps to plot in a different colour
##' @param late.col colour in which to plot final `late.num` time steps
##' @param early.col colour in which to plot earlier time step points
##' @param early.col.lines colour in which to plot earlier time step points
##' @param start first time step (must be 1)
##' @param pt.type `type` value for `points()`
##' @param par.mgp `par("mgp")` values
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   plot_time_series(NY_lags_example$N_t, X.or.N = "N")
##'   plot_time_series(NY_lags_example$Y_t, X.or.N = "X")
##' }
plot_time_series <- function(values,
                             X.or.N,
                             par.mar.ts = c(3, 3, 1, 1),
                             max_time = NULL,
                             t.axis.range = NULL,
                             last.time.to.plot = NULL,
                             late.num = 3,
                             late.col = "red",
                             early.col = "black",
                             early.col.lines = "lightgrey",
                             start = 1, # may not work for others
                             pt.type = "p",
                             par.mgp = c(1.5, 0.5, 0)
                             ){
  stopifnot(start == 1)

  if(is.null(max_time)) max_time <- length(values)

  if(is.null(last.time.to.plot)){
    if(X.or.N == "N"){
       last.time.to.plot <- max(which(!is.na(values)))} else {
                                                      last.time.to.plot <-
                                                        length(values)
                                                      # Want the NA included for
                                                      #  X, but N has an extra
                                                      #  time step added
                                                      }
  }


  if(is.null(t.axis.range)) {
      t.axis.range <- c(start, max_time)
  }

  par("mgp" = par.mgp) # first val sets axis title dist (mex units)
                       #  to axes, second sets labels
  par(pty = "m")                 # maximal plotting region, not square
                                 #  like for phase plots
  par(mar = par.mar.ts)

  col.plot = c(rep(early.col,
                   max(c(0, last.time.to.plot - late.num))),
               rep(late.col,
                   min(c(last.time.to.plot, late.num))) )   # colours of points
  col.plot.lines = col.plot                            # colours of lines
  col.plot.lines[col.plot.lines == early.col] = early.col.lines
  pch.plot = (col.plot == early.col) * 1 + (col.plot == late.col) * 16
                                        # filled circles for latest
  pch.plot[length(pch.plot)] = 8     # latest one a star


  if(X.or.N == "N"){
    Nt.max.abs = max( abs( range(values[start:max_time],
                                 na.rm=TRUE) ) )
    Nt.axes.range = c(0, Nt.max.abs*1.04)    # Expand else points can hit edge

    plot(0, 0,
         xlab = expression("Time, t"),
         ylab = expression("N"[t]),
         xlim = c(0, max(t.axis.range)),
         ylim = Nt.axes.range,
         type = "n",                           # empty plot
         main = paste0("Time t=", last.time.to.plot))
  } else {
    XtLoc = -0.05 * max(t.axis.range)  # location to plot Xt on a vertical line,
    Xt.max.abs = max(abs( range(values[start:max_time],
                                na.rm=TRUE) ),
                     abs( range(values[start:max_time],
                                na.rm=TRUE)) )
    Xt.axes.range = c(-Xt.max.abs, Xt.max.abs)

    plot(0, 0,
         xlab = expression("Time, t"),
         ylab = expression("Y"[t]),
         xlim = c(XtLoc, max(t.axis.range)),
         ylim = Xt.axes.range,
         type = "n")                           # empty plot
    abline(v = 0.5*XtLoc, col="black")
  }

  iii = last.time.to.plot             # use iii since simpler
  if(iii > 1.5){
    segments(start:(iii-1),
             values[start:(iii-1)],    # dplyr::pull(Nx.lags.use[start:(iii-1), "Y_t"]),
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
  invisible()
}

##' Plot 2d phase plot of `N_t` vs `N_{t-1}` or `Y_t` vs `Y_{t-1}`
##'
##' Shows the values in two-dimensional phase space (first differenced or not)
##'   with lag of 1. Cobwebbing shows how one point iterates to the next point
##'   in the time series.
##'
##' @param values vector of `N_t` or `Y_t` values
##' @param X.or.N "N" if raw non-differenced data, "X" for differenced data
##' @param par.mar.phase `par(mar)` values
##' @param axis.range range of axes, if NA then calculated from values
##' @param start first time step to plot (currently must be 1)
##' @param last.time.to.plot last time value of N[t] to use when plotting, so
##'   final Y[t] used will be Y[t-1] (since Y[t] uses N[t+1])
##' @param pt.type `type` value for `points()`
##' @param cobwebbing if TRUE then add cobwebbing lines to phase plot
##' @param late.col colour in which to plot final `late.num` time steps
##' @param early.col colour in which to plot earlier time step points
##' @param early.col.lines colour in which to plot earlier time step points
##' @param late.num final number of `N[t]` time steps to plot in a different colour
##' @return Plots figure to current device
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   plot_phase_2d(NY_lags_example$N_t, X.or.N = "N")
##'   plot_phase_2d(NY_lags_example$Y_t, X.or.N = "X")
##' }
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
                          late.num = 3
                          ){

  if(is.na(axis.range)) {
    axis.range <- c(min(0, min(values, na.rm = TRUE)),
                    max(values, na.rm = TRUE))   # for N_t or Y_t
  }
  if(is.null(last.time.to.plot)) last.time.to.plot <-
                                   max(which(!is.na(values)))

  stopifnot(start == 1)
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
         xlab = expression("Y"[t-1]),
         ylab = expression("Y"[t]),
         xlim = axis.range,
         ylim = axis.range,
         type = "n")
    values.to.plot <- values[start:(last.time.to.plot-1)] # Not use Y[last..]
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
        # dplyr::pull(Nx.lags.use[start:(iii-1), "N_tmin1"]),
        # dplyr::pull(Nx.lags.use[start:(iii-1), "N_t"]),
        # dplyr::pull(Nx.lags.use[(start+1):iii, "N_tmin1"]),
        # dplyr::pull(Nx.lags.use[(start+1):iii, "N_t"]),
        # not fully tested:
        pbsLag(values.to.plot)[-length(values.to.plot)],   # N_{t-1}
        values.to.plot[-length(values.to.plot)],
        pbsLag(values.to.plot)[-1],
        values.to.plot[-1],
        col = col.plot.lines) # lines() will not use vector of col
    }
  }
  if(last.time.to.plot > 1.5){
    points(pbsLag(values.to.plot),
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
#               xvals = rep( dplyr::pull(Nx.lags.use[start:iii, "Y_t"]), each = 2)
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
#               segments(dplyr::pull(Nx.lags.use[start:(iii-1), "Y_tmin1"]),
#                        dplyr::pull(Nx.lags.use[start:(iii-1), "Y_t"]),
#                        dplyr::pull(Nx.lags.use[(start+1):iii, "Y_tmin1"]),
#                        dplyr::pull(Nx.lags.use[(start+1):iii, "Y_t"]),
#                        col = col.plot.lines)
#            }
#        }
#      if(iii > 1.5)
#        {
#          points(dplyr::pull(Nx.lags.use[start:iii, "Y_tmin1"]),
#                 dplyr::pull(Nx.lags.use[start:iii, "Y_t"]),
#                 type = pt.type, pch = pch.plot,
#                 col = col.plot)           # start row has NA's, get ignored
#        }
      # legend("topleft", legend=paste("Time", iii), border = NULL)
  invisible()
}

##' Plot 3d phase plot of `Y_t` vs `Y_{t-1}` vs `Y_{t-2}`
##'
##' Shows the first-differenced values in three-dimensional phase space, as
##'   points showing lags of 0, 1 and 2. Cobwebbing shows how one point iterates
##'   to the next point in the time series.
##'
##' @param obj `pbsEDM` object
##' @param par.mgp.3d `par(mgp)` values
##' @param par.mai.3d `par(mai)` values
##' @param par.mar.3d `par(mar)` values
##' @param x.lab x axis label
##' @param y.lab y axis label
##' @param z.lab z axis label
##' @param last.time.to.plot last time value of N[t] to use when plotting, so
##'   final Y[t] used will be Y[t-1] (since Y[t] uses N[t+1])
##' @param axis.range range of axes, if NA then calculated from values; all
##'   three axes have the same range
##' @param late.num final number of `N_t` time steps to plot in a different
##'   colour; not showing `N_t` here but keeping consistency with other plots
##' @param pt.type `type` value for `points()`
##' @param late.col colour in which to plot final `late.num` time steps
##' @param early.col colour in which to plot earlier time step points
##' @param early.col.lines colour in which to plot earlier time step points
##' @param axes.col colour of orthogonal origin axes that go through (0, 0, 0)
##' @param par.mar.phase `par(mar)` to reset to after plotting 3d figure
##' @param par.mgp `par(mgp)` to reset to after plotting 3d figure
##'
##' @return Plots figure to current device
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'  aa <- pbsEDM(NY_lags_example,
##'               lags = list(N_t = 0:1),
##'               first_difference = TRUE)
##'   plot_phase_3d(aa)
##' }
plot_phase_3d <- function(obj,
                          par.mgp.3d = c(3, 10, 0),
                          par.mai.3d = c(0.1, 0.1, 0.1, 0.1),
                          par.mar.3d = c(3, 0, 0, 0),
                          x.lab = expression("Y"[t-2]),
                          y.lab = expression("Y"[t-1]),
                          z.lab = expression("Y"[t]),
                          last.time.to.plot = NULL,
                          axis.range = NA,
                          late.num = 3,
                          pt.type = "p",
                          late.col = "red",
                          early.col = "black",
                          early.col.lines = "lightgrey",
                          axes.col = "darkblue",
                          par.mar.phase = c(3, 0, 1, 0),  # to reset for normal figs
                          par.mgp = c(1.5, 0.5, 0)
                          ){
  if(is.na(axis.range)) {
    axis.range <- c(min(0, min(obj$X_observed, na.rm = TRUE)),
                    max(obj$X_observed, na.rm = TRUE))
  }

  if(is.null(last.time.to.plot)) last.time.to.plot <- length(obj$X_observed)

  start = 1     # currently only works for 1

  if(last.time.to.plot > 2){
    values.to.plot <- obj$X_observed[start:(last.time.to.plot-1)]   # can't use N[last...])
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
#                                                            "Y_tmin2"]),
#                                    dplyr::pull(Nx.lags.use[start:iii,
#                                                            "Y_tmin1"]),
  #                                    dplyr::pull(Nx.lags.use[start:iii, "Y_t"]) )
  # maybe this is okay with NA's:
#  proj.pts = scat$xyz.convert(pbsLag(obj$X_observed,
#                                     2)[start:iii],  # "Y_tmin2", index wrong?
#                              pbsLag(obj$X_observed,
#                                     1)[start:iii],  # "Y_tmin1"
  #                              pbsLag(obj$X_observed)[start:iii])  # "Y_t"
  if(all(!is.na(values.to.plot))){
    proj.pts = scat$xyz.convert(pbsLag(values.to.plot,
                                       2),  # "Y_tmin2"
                                pbsLag(values.to.plot,
                                       1),  # "Y_tmin1"
                                values.to.plot)  # "Y_t"
    if(last.time.to.plot > 3.5){
      segments(proj.pts$x[1:(last.time.to.plot - start)],
               proj.pts$y[1:(last.time.to.plot - start)],
               proj.pts$x[2:(last.time.to.plot - start + 1)],
               proj.pts$y[2:(last.time.to.plot - start + 1)],
               col = col.plot.lines) # lines() will not use vector
    }

    # The points
    if(last.time.to.plot > 2.5){
      scat$points3d(pbsLag(values.to.plot,
                           2),  # "Y_tmin2"
                    pbsLag(values.to.plot,
                           1),  # "Y_tmin1"
                    values.to.plot,  # "Y_t"
                    type = pt.type,
                    pch = pch.plot,
                    col = col.plot)
    }
  }

  par(mar = par.mar.phase)   # scatterplot3d changes mar
  par(mgp = par.mgp)         #  so set back to usual for 2d figures
  invisible()                   # else returns the above line
}

##' 2-d phase plot of x_t v x_{t-1} with coloured points to explain EDM
##'
##' Highlights a point to be projected, its nearest `E+1` neighbours, and then
##' draw arrows to show where they go and so where the projection goes. Very
##' useful for understanding and checking what EDM is doing. Call this multiple
##' times using TODO to make a movie. Type `plot_explain_edm_movie` for example calls.
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param tstar the time index, `t*` for `x(t*)` of the target point for which to make a
##'   projection from
##' @param x.lab label for x axis
##' @param y.lab label for x axis
##' @param main main title
##' @param tstar.col colour for `x(t*)`
##' @param tstar.pch pch value (point style) for `x(t*)`
##' @param tstar.cex cex value (size) for `x(t*)`
##' @param neigh.plot if TRUE then highlight neighbours to `x(t*)`
##' @param neigh.proj if TRUE then highlight projections of neighbours to `x(t*)`
##' @param pred.plot if TRUE then highlight forecasted value `x(t*+1)`
##' @param pred.rEDM if TRUE then show predictions from `rEDM` for
##'   `NY_lags_example` saved values; obj must be `NY_lags_example` TODO do a test as
##'   may not work
##' @param true.val if TRUE then plot the true value of `x(t*+1)`
##' @param legend.plot if TRUE then do a legend and print value of `t*`
##' @return single plot that explains one part of EDM, link together in a movie
##'   using `pbs_explain_edm_movie()`
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM(NY_lags_example,
##'               lags = list(N_t = 0:1),
##'               first_difference = TRUE)
##'   plot_explain_edm(aa,
##'                    tstar = 15,
##'                    main = paste0(
##'                            "See where neighbours go"),
##'                    tstar.pch = 19,
##'                    tstar.cex = 1.2,
##'                    neigh.plot = TRUE,
##'                    neigh.proj = TRUE)
##' # Type plot_explain_edm_movie to see other examples
##' }
plot_explain_edm <- function(obj,
                             tstar,
                             x.lab = expression("Y"[t-1]),
                             y.lab = expression("Y"[t]),
                             main = "All the points in lagged space",
                             tstar.col = "blue",
                             tstar.pch = 1,
                             tstar.cex = 1,
                             neigh.plot = FALSE,
                             neigh.proj = FALSE,
                             pred.plot = FALSE,
                             pred.rEDM = FALSE,
                             true.val = FALSE,
                             legend.plot = TRUE){
  if(pred.rEDM){
    testthat::expect_equal(obj$X_observed, NY_lags_example$Y_t)
    #  ideally want this message:
    #  stop("pred.rEDM can only be TRUE when plotting NY_lags_example")
  }

  par(pty="s")
                             # TODO maybe. Andy NOT CHANGED THESE (yet) in notation update
  Xt <- obj$X[, "N_t_0"]      # Should change Nt in pbsEDM() as confusing. Luke
                              # was still tweaking that. Think here want Xt ->
                              # Y_t, and Xtmin1 -> X_tmin1
  Xtmin1 <- obj$X[, "N_t_1"]

  Xt.max.abs = max(abs( range( c(Xt,
                                 Xtmin1,
                                 obj$X_forecast),
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

  psi_vec <- obj$neighbour_index[tstar,]     #  vector of indices of points closest
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
               #         sh.lwd=2, width = 1, sh.col="red")
    }
  }

  # plot forecast value
  if(pred.plot){
    points(Xt[tstar],
           obj$X_forecast[tstar + 1],    # CHECK not tstar+1, think it's correct
           col = tstar.col,
           pch = 8)

    igraph:::igraph.Arrows(Xtmin1[tstar],
                           Xt[tstar],
                           Xt[tstar],
                           obj$X_forecast[tstar + 1],
                           curve = 1,
                           size = 0.7,
                           h.lwd = 2,
                           sh.lwd = 2,
                           width = 1,
                           sh.col = "darkgrey")

    abline(h = Xt[psi_vec + 1],
           col = "lightgrey")
  }

  # show predicted value from rEDM
  if(pred.rEDM){
    points(Xt[tstar],
           dplyr::pull(NY_lags_example[tstar+1, "rEDM.pred"]),
           col = "red",
           cex = 1.5,
           lwd = 2)
  }

  # plot the true value of xt(t^*+1)
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
##' Use with gifski in .Rmd - see `analyse_simple_time_series` vignette.
##'
##' @param obj list object of class `pbsEDM`, an output from `pbsEDM()`
##' @param ... extra values to use in calls to `plot_explain_edm()`; needs to include tstar
##' @return movie (if used with gifski) to explain EDM; if combined with
##'   `par(mfrow=c(4,2))` (and not gifski) then will give multiple panels, but
##'   borders etc have not been set up properly for this and will need adjusting
##' @export
##' @author Andrew Edwards
##' @examples
##' \donttest{
##'   aa <- pbsEDM(NY_lags_example,
##'               lags = list(N_t = 0:1),
##'               first_difference = TRUE)
##'   plot_explain_edm_movie(aa, tstar = 15) # will only show last panel; see
##'                                          # `analyse_simple_time_series`
##'                                          #  vignette.
##' }
plot_explain_edm_movie <- function(obj,
                                   ...){
  # Need to specify tstar, but plot_explain_edm() will give automatic error anyway
  # if(!exists("tstar")){
  #   stop("Need to specify tstar")
  # }

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
                     "and take weighted average of Y_t to be predicted value"),
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

#' Plot Method for `pbsSim`
#'
#' @param x [pbsSim()]
#' @param ... Other arguments.
#'
#' @return Panel plot.
#' @export
#'
plot.pbsSim <- function (x, ...) {
  par(mfrow = c(3, 1), mar = c(4, 4, 1, 0), oma = c(2, 1, 1, 1))
  plot(x[, "producers"], type = "l", ylab = "Producers")
  plot(x[, "prey"], type = "l", ylab = "Prey")
  plot(
    x[, "predators"],
    type = "l",
    xlab = "Time step",
    ylab = "Predators")
  par(mfrow = c(1, 1))

}

##' Plot the size of the library for a given `T` as a function of `E` and `tstar`
##'
##' Show a coloured contour plot of how the library size changes for a
##' univariate time series with sample size `T` as chosen embedding dimension
##' `E` and target index `tstar` vary.
##'
##' @param T integer of the sample size of a univariate time series (time series
##'   itself not needed)
##' @param E_vec vector of integers to show for the embedding dimension
##' @param tstar_vec vector of integers to show for the focal point time index
##'   `tstar`, `t*` in the manuscript, from which projections are made from
##' @param annotate logical whether to add numbers along the middle of the plot
##'   (avoids the need for a colorbar, which was fiddly to do, see attempts
##'   deleted in 8f567ff, while also having a grey background), and add text to
##'   explain the grey areas.
##' @param annotate_cex text size for main annotation
##' @param annotate_tstar tstar value at which to add the annotated numbers; if
##'   NULL then is `max(tstar_vec)/2`.
##' @param annotate_extra_cex text size for extra smaller numbers in boxes at
##'   the top
##' @param words_pos vector of positioning for words in grey areas, given by `E`
##'   and `tstar` position for `Early values' text (text is placed to the right
##'   of these), then for `Late values' text (text is centred around these.
##' @param words_cex text size for above words
##' @param mgp_vals mgp values to use, as in `plot(..., mgp = mgp_vals)`.
##' @return plots the contour plot and returns the matrix (`C` in manuscript) of calculated library sizes
##' @export
##' @author Andrew Edwards
##' @examples
##' \dontrun{
##' # For manusript doing:
##' postscript("library_size.eps",
##'             height = 6,
##'             width = 6,
##'             horizontal=FALSE,
##'             paper="special")
##' plot_library_size()
##' dev.off()
##' }
plot_library_size <- function(T = 50,
                              E_vec = 2:10,
                              tstar_vec = 1:50,
                              annotate = TRUE,
                              annotate_cex = 1,
                              annotate_tstar = 23,
                              annotate_extra_cex = 0.7,
                              words_pos = c(6, 3.5, 6, 49.6),
                              words_cex = 0.9,
                              mgp_vals = c(2, 0.5, 0)){
  stopifnot(length(T) == 1,
            length(E_vec) > 1,
            min(E_vec) > 1,
            length(tstar_vec) > 1)

  stopifnot("T is too small relative to max(E_vec); change code if you want to try exceptions" =
              T - 2 * (max(E_vec) + 1) >= 0)

  if(is.null(annotate_tstar)){
    annotate_tstar <- max(tstar_vec)/2
  }

  C <- matrix(NA,
              nrow = length(E_vec),
              ncol = length(tstar_vec))

  for(i in 1:length(E_vec)){
    E <- E_vec[i]
    if(T - E - 2 >= E){
      for(tstar in E:(T - E - 2)){
        j <- which(tstar_vec == tstar)
        C[i, j] <- T - 2 * (E + 1)
      }
    }
    for(tstar in (T - E - 1):(T - 2)){    # E > 1 so always valid
      j <- which(tstar_vec == tstar)
      C[i, j] <- tstar - E
    }
    # And tstar <= E-1 and tstar >= T-1 remain as NA
  }

  image(E_vec,
        tstar_vec,
        C,
        xlim = range(E_vec) + c(-0.5, 0.5),   # colours are correctly around integers
        ylim = rev(range(tstar_vec) + c(-0.5, 0.5)),
        xaxs = "i",
        yaxs = "i",         # sets exact axis range
        xlab = expression(paste("Embedding dimension, ", italic(E))),
        ylab = expression(paste("Focal time, ", italic(t), "*")),
        mgp = mgp_vals)

  rect(0,
       0,
       max(E_vec) + 2,
       max(tstar_vec) + 2,
       col = "darkgrey")    # So NA's come out grey

  image(E_vec,
        tstar_vec,
        C,
        add = TRUE,
        col = rev(hcl.colors(max(C, na.rm=TRUE) - min(C, na.rm=TRUE) + 1,
                         "Spectral")))

  # Add extra unlaballed tickmarks
  box()
  axis(1, E_vec, tcl = -0.4, labels = rep("", length(E_vec)))
  axis(2, tstar_vec, tcl = -0.2, labels = rep("", length(tstar_vec)))
  every_five <- tstar_vec[which(tstar_vec %% 5 == 0)]
  axis(2, every_five, tcl = -0.4, labels = rep("", length(every_five)))

  if(annotate){
    annotate_main <- C[1:length(E_vec),
                       which(tstar_vec == annotate_tstar)]  # main annotation
    text(E_vec,
         annotate_tstar,
         annotate_main,
         cex = annotate_cex)

    annotate_extra <- which(C > max(annotate_main),
                            arr.ind = TRUE)          # extra small annotations

    for(k in 1:nrow(annotate_extra)){
      text(E_vec[annotate_extra[k, "row"]],
           tstar_vec[annotate_extra[k, "col"]],
           C[annotate_extra[k, "row"],
             annotate_extra[k, "col"]],
           cex = annotate_extra_cex)
    }

    text(words_pos[1],
         words_pos[2],
         expression(paste("Early values of ",
                          italic(t),
                          "* not possible")),
         pos = 4,
         cex = words_cex)
    text(words_pos[3],
         words_pos[4],
         expression(paste("Late values of ",
                          italic(t),
                          "* not possible")),
         cex = words_cex)
  }

  return(C)
}
