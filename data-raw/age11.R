# Results from 2026 SalmonPrize competition for Columbia Bonneville,
# to use as a small and large example output from MVE to create a summary figure for the lags.

# On Andy's computer:

# Only has 11 lag combinations in the top subsets
age11_res <- readRDS(paste0(here::here(),
                            "/../sockeyePrize/report-2026/columbia/bonneville/age11_res.rds"))

# 44 lag combinations in the top subsets
age12_res <- readRDS(paste0(here::here(),
                            "/../sockeyePrize/report-2026/columbia/bonneville/age12_res.rds"))

# Write to data/
usethis::use_data(age11_res, overwrite = TRUE)
usethis::use_data(age12_res, overwrite = TRUE)
