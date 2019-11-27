#' Create a matrix of lags of an original time series vector
#'
#' @param time_series A numeric vector with possible NA values
#' @param embed_dim An integer embedding dimension
#' @param lag_size The number of time steps separating successive lags
#'
#' @return matrix of successively lagged columns
#' @export
#'
#' @examples
pbs_make_lags <- function(time_series, embed_dim, lag_size) {
  ts_index <- tibble::tibble(index = seq_len(NROW(time_series)))
  lapply(
    X = seq_len(embed_dim) - 1,
    FUN = function(X, ts, n) dplyr::lag(ts, X * n),
    ts = time_series,
    n = lag_size
  ) %>% 
    dplyr::bind_cols(ts_index) %>%
    drop_na() %>%
    right_join(ts_index, by = c("index")) %>%
    select(-index) %>%
    as.matrix()
}