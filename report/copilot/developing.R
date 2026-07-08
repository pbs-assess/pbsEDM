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

# Great, the function works as specified. Now add an extra argument, `order`, to specify
# the order of the new strings in the output. For example, `order = c("age11",
# "pink")` would ensure that the strings are returned in the order
# `age11_lag_1`, `age11_lag_2, `age11_lag_3`, `pink_lag_1`. For remaining
# variable names that are not specified in `order`, just put them in
# alphabetical order. If `order` is `NULL`, then return in alphabetical
# order. The `lag_1` etc. should always be in numerical order of the number.

# Please suggest a plan.
