
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pbsEDM

<!-- badges: start -->

[![R-CMD-check](https://github.com/pbs-assess/pbsEDM/workflows/R-CMD-check/badge.svg)](https://github.com/pbs-assess/pbsEDM/actions)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Coverage
status](https://codecov.io/gh/pbs-assess/pbsEDM/branch/master/graph/badge.svg)](https://codecov.io/github/pbs-assess/pbsEDM?branch=master)
<!-- badges: end -->

An R package to implement the delay embedding methods of Empirical
Dynamic Modeling

## Installation

``` r
devtools::install_github("pbs-assess/pbsEDM")
```

## To reproduce and save the figures for the manuscript submitted to Methods in Ecology and Evolution

```E_results <- pbsEDM_Evec(NY_lags_example$N_t)```

Figure 1 - ```plot_pbsEDM_Evec_save(E_results)```

Figure 2 -  ```plot_explain_edm_save(E_results[[1]])```

Figure 3 - ```plot_rho_Evec_save()```

Figure 4 - ```plot_library_size_save()```

Animated figures for the Appendix:

Figure S.1 - 

```plot_pbsEDM_Evec_movie_save(E_results)```

Figure S.2 - 


```plot_explain_edm_movie_save(E_results[[1]])```

Figure S.3 - 

```plot_explain_edm_all_tstar_movie_save(E_results[[1]])```

## Workflow During Development

Both pushing to same repo, so always check on GitHub first, then (if
necessary)

    git fetch
    git rebase

If get conflicts then see some notes at
<https://github.com/pacific-hake/hake-assessment> (for which we are also
pushing to a single repository), scroll down to the GitHub workflow
section.
