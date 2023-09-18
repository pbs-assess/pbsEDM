#' Weight Single View Embeddings
#'
#' @param forecasts [list()]
#' @param metric [character()]
#' @param weight [character()]
#' @param n_weight [numeric()]
#'
#' @importFrom rlang .data
#' @importFrom rlang :=
#'
#' @return [list()]
#' @export
#'
weight_single_view_embeddings <- function (forecasts,
                                           metric,
                                           weight,
                                           n_weight) {
  # Define ranks
  ranks <- forecasts %>%
    dplyr::bind_rows(.id = "ssr") %>%
    dplyr::mutate(ssr = as.numeric(.data$ssr)) %>%
    dplyr::arrange(.data$time, .data[[metric]]) %>%
    dplyr::group_by(.data$time) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$ssr, .data$time) %>%
    dplyr::group_by(.data$ssr) %>%
    dplyr::mutate(lag_rank = dplyr::lag(.data$rank, n = 1L)) %>%
    dplyr::ungroup()
  # Define forecast
  forecast <- ranks %>%
    dplyr::filter(.data$lag_rank <= n_weight) %>%
    dplyr::arrange(.data$time, .data$lag_rank) %>%
    dplyr::group_by(.data$time) %>%
    dplyr::mutate(forecast = mean(.data$forecast, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.data$time, .keep_all = TRUE) %>%
    dplyr::select(.data$set, .data$time, .data$observed,.data$forecast)
  # Define summary
  summary <- ranks %>%
    dplyr::filter(.data$lag_rank <= n_weight) %>%
    dplyr::arrange(.data$time, .data$lag_rank) %>%
    dplyr::filter(.data$set == 1)
  # Define hindsight
  hindsight_ssr <- ranks %>%
    dplyr::filter(.data$time == max(.data$time, na.omit = TRUE)) %>%
    dplyr::filter(.data$rank == 1) %>%
    dplyr::pull(.data$ssr)
  hindsight <- ranks %>%
    dplyr::filter(.data$ssr == hindsight_ssr) %>%
    dplyr::filter(.data$set == 1)
  # Define observed-forecast matrix
  m <- matrix(c(forecast$observed, forecast$forecast), ncol = 2L)
  # Define results
  results <- forecast %>%
    dplyr::mutate(
             mre = runner::runner(m, f = matric, fun = mre),
             !!metric := runner::runner(m, f = matric, fun = get(metric))
           ) %>%
    dplyr::filter(.data$set == 1) %>%
    dplyr::select(.data$time:.data[[metric]])
  # Return
  return(
    list(
      ranks = ranks,
      summary = summary,
      hindsight = hindsight,
      results = results
    )
  )
}
