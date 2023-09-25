#' Multiview Embedding (Luke's original function).
#'
#' See Andy's `multiview_embedding()` for what he's using now.
#'
#' From Luke's email (TODO tidy up when done):
#' Good to hear from you! Yes, what you describe is the main idea behind
#' multiview embedding (MVE). It uses the `2^N - 1` subsets of the superset
#' state space reconstruction of covariates and covariate lags, giving all
#' possible embedding dimensions and embedding combinations that you're looking
#' for. [Andy: Ye and Sugihara 2016 exclude those that have no 0 lag variables,
#' since presumably they're just shifted - keep in for now and check results
#' agree maybe].
#' Andy asked: is single-view embedding just standard EDM, and Luke said:
#' I think MVE was originally a generalization of S-mapping [Andy: no, for MVE
#' they just use the single nearest neighbour in each view], but we've
#' implemented it as a generalization of EDM. So in our case, single-view
#' embedding would be EDM, but other authors might hear 'single-view embedding'
#' and think 'S-mapping'.
#' Worth double-checking because this is just from memory ;)
#'  Andy: TODO don't think that's quite right, checking
#' the code.
#' TODO TODO TODO check the definitions; stick with mve for now to get it all
#' working.
#'
#' @param data [matrix()] or [data.frame()] with named [numeric()] columns
#' @param response [character()] column name of the response variable in
#'   \code{data}
#' @param lags [list()] of a named vector of lags for each explanatory
#'   variable.
#' @param index [integer()]
#' @param buffer [integer()] number of forecasts prior to \code{index}
#' @param window [integer()] forecast metric moving window width
#' @param metric [character()]
#' @param beyond [logical()]
#' @param weight TBD
#' @param n_weight [integer()]
#' @param cores [integer()] number of cores for parallel computation.
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
multiview_embedding_original <- function (data,
                 response,
                 lags,
                 index = 50L,
                 buffer = 10L,
                 window = integer(0),
                 metric = "rmse",
                 beyond = FALSE,
                 weight = NULL,
                 n_weight = ceiling(sqrt(2^length(unlist(lags)))),
                 cores = NULL) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(index, lower = 20, upper = nrow(data), len = 1)
  checkmate::assert_integerish(buffer, lower = 1, upper = 10, len = 1)

  # Create subset lags ---------------------------------------------------------

  subset_lags <- create_subset_lags(lags)

  # Apply running empirical dynamic modeling -----------------------------------

  if (is.null(cores)) {
    # Apply in sequence
    forecasts <- lapply(
      subset_lags,
      FUN = single_view_embedding,
      data = data,
      response = response,
      index = index,
      buffer = buffer,
      window = window,
      metric = metric,
      beyond = beyond,
      superset = lags
    )
  } else {
    if (.Platform$OS.type == "unix") {
      # Apply in parallel via forking
      forecasts <- parallel::mclapply(
                               subset_lags,
                               FUN = single_view_embedding,
                               data = data,
                               response = response,
                               index = index,
                               buffer = buffer,
                               window = window,
                               metric = metric,
                               beyond = beyond,
                               superset = lags,
                               mc.cores = cores
                             )
    } else {
      # Apply in parallel via sockets
      # cl <- parallel::makeCluster(cores)
      # parallel::clusterEvalQ(cl, library()) # TODO: packages needed
      # results_list <- parallel::parLapply()
      stop("eedm parallel via sockets not yet implemented R/functions.R")
    }
  }

  # Weight redm forecasts ------------------------------------------------------

  weighted <- weight_single_view_embeddings(forecasts, metric, weight, n_weight)

  # Return results object ------------------------------------------------------

  return(
    structure(
      list(
        data = data,
        observed = c(dplyr::pull(data, response), NA),
        forecast = c(rep(NA_real_, index - 1), weighted$results$forecast),
        response = response,
        lags = lags,
        index = index,
        buffer = buffer,
        window = window,
        metric = metric,
        beyond = beyond,
        n_weight = n_weight,
        raw_forecasts = forecasts,
        ranks = weighted$ranks,
        summary = weighted$summary,
        hindsight = weighted$hindsight,
        results = weighted$results
      ),
      class = "multiview_embedding"
    )
  )
}
