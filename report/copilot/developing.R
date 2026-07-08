# Instructions for copilot are given as comments (i.e. after #), with R code to
# run not commented.

# This is an R package, for which I want to add a new function, as described in
# the comments below.

load_all()

lags_example <- age11_res$subset_lags[1:10]

# lags_example is a list, where each element contains one or more variable names
# and one or more numbers (know as lags) associated with that variable name.
# Create a function called `create_named_lags()` in R/create-named-lags.R. It will input the list object and return a vector of strings denoting each unique combination of name and
# lag.

# For example, the first element in the list is `age11` with a lag of 1, and so the result
# for that element should be `age11_lag_1`.
# The third element is `age11` with
# lags of 1 and 2, and so the result for that element should be two values:
# `age11_lag_1` and `age11_lag_2`.

# Please suggest a plan.
