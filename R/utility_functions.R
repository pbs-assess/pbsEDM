#' Create a Tibble of Lagged Copies of a Named Tibble Column
#'
#' @param tbl A tibble with a named column to be lagged
#' @param col_name The name of a column in tbl
#' @param embed_dim An integer embedding dimension corresponding to the number
#'     of lagged columns to create
#' @param lag_size The number of time steps separating successive lags
#'
#' @return A tibble of successively lagged columns
#' 
#' @importFrom magrittr %>%
#' @export
#'
#' @examples pbs_make_lags(tibble::tibble(x = 1:15), "x", 5, 1)
#' 
pbs_make_lags <- function(tbl,
                          col_name,
                          embed_dim,
                          lag_size) {
  # Check arguments
  stopifnot(
    tibble::is_tibble(tbl),
    col_name %in% names(tbl),
    is.numeric(embed_dim),
    round(embed_dim) == embed_dim,
    round(embed_dim) >= 1L,
    is.numeric(lag_size),
    round(lag_size) == lag_size,
    round(lag_size) >= 1L
  )
  
  # Create time_series vector and ts_index tibble
  time_series <- dplyr::pull(tbl, col_name)
  ts_index <- tibble::tibble(index = seq_along(time_series))
  
  # Return lagged tibble
  lapply(
    X = seq_len(embed_dim) - 1,
    FUN = function(X, ts, n) dplyr::lag(ts, X * n),
    ts = time_series,
    n = lag_size
  ) %>%
    dplyr::bind_cols(ts_index) %>%
    tidyr::drop_na() %>%
    dplyr::right_join(ts_index,
                      by = c("index")) %>%
    dplyr::select(-index) %>%
    magrittr::set_colnames(paste0("lag_", as.character(seq_len(embed_dim) - 1)))
}
