# As for sizeSpectra package, need these to avoid warnings related to dplyr
# commands (e.g. referring to the column names within dplyr::filter() ).
# Copied from the warning given by check() (that puts them alphabetical, then
# maybe add some more at the end).
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
           "Nx_lags_orig",
           "Xt",
           "Xtmin1",
           "d",
           "my.pred",
           "rEDM.pred"
  ))
}
