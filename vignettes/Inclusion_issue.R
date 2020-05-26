## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6
)

## ----setup--------------------------------------------------------------------
library(pbsEDM)

## ----fullresults--------------------------------------------------------------
as.data.frame(Nx_lags_orig)

## -----------------------------------------------------------------------------
Nx_lags_orig[100, c("rEDM.pred", "my.pred")]

## -----------------------------------------------------------------------------
plotPanelMovie.df2(only.final.plot=TRUE,
                   open.pdf=FALSE)

## -----------------------------------------------------------------------------
eps = 0.00001      # How different the predictions can be
plot(Nx_lags_orig$rEDM.pred,
     Nx_lags_orig$my.pred,
     xlab = "rEDM predictions",
     ylab = "Andy predictions")
abline(a=0,
       b=1,
       col="grey")
# Colour in red the ones more than eps away
different <- dplyr::filter(Nx_lags_orig,
              abs(pred.diff) > eps)
points(different$rEDM.pred,
       different$my.pred,
       col = "red",
       pch = 20)


## -----------------------------------------------------------------------------
different

## ----pbsEDMcalc---------------------------------------------------------------
# pbs_calc <- pbs_edm(Nx_lags_orig,
#                     lags = list(Xt = c(0:1))) # A tibble (contains lists)
# 
# testthat::expect_equal(Nx_lags_orig$Xt,
#                        pbs_calc$observations$observations) # check the indexing is the same
# 
# pbs_pred <- pbs_calc$forecasts$forecasts
# 
# plot(Nx_lags_orig$rEDM.pred,
#      pbs_pred,
#      xlab = "rEDM predictions",
#      ylab = "pbsEDM predictions")
# abline(a=0,
#        b=1,
#        col="grey")

# Colour in red the ones more than eps away
# different <- dplyr::filter(Nx_lags_orig,
#              abs(pred.diff) > eps)
#points(different$rEDM.pred,
#       different$my.pred,
#       col = "red",
#       pch = 20)

