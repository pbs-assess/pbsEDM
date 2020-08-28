
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pbsEDM

[![Travis build
status](https://travis-ci.org/luke-a-rogers/pbsEDM.svg?branch=master)](https://travis-ci.org/luke-a-rogers/pbsEDM)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Coverage
status](https://codecov.io/gh/luke-a-rogers/pbsEDM/branch/master/graph/badge.svg)](https://codecov.io/github/luke-a-rogers/pbsEDM?branch=master)

An R package to implement the delay embedding methods of Empirical
Dynamic Modeling

## Installation

``` r
devtools::install_github("luke-a-rogers/pbsEDM")
```

## Notation

Where possible we have followed the useful notation in Deyle et al. (2013), namely:

| $t$ | Time, taking integer values 1, 2, 3, ... |
| $N(t)$ | Population at time $t$ |
| $X(t)| First difference: $X(t) = N(t+1) - N(t)$, $t = 1, 2, 3, ...$ |


## Functionality

TBA

## Examples

TBA

## Workflow During Development

Both pushing to same repo, so always check on GitHub first, then (if
necessary)

    git fetch
    git rebase

If get conflicts then see some notes at
<https://github.com/pacific-hake/hake-assessment> (for which we are also
pushing to a single repository), scroll down to the GitHub workflow
section.
