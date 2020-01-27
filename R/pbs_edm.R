#' Use Empirical Dynamic Modeling to Perform Out-of-Sample Forecasting 
#'
#' @param data_frame A data frame with named columns (list of numeric vectors)
#' @param lags Named list of numeric vectors giving lags to use for each
#'     column. List names must match column names. (list of numeric vectors)
#' @param forecast_distance Forecast distance (numeric scalar)
#' @param symmetric_exclusion Symmetric exclusion radius? (logical scalar)
#' @param include_stats Return forecast stats? (logical scalar)
#' @param include_forecasts Return forecast and observations? (logical scalar)
#' @param include_neighbours Return neighbours? (logical scalar)
#'
#' @return A list with results, forecasts, neighbours, distances, and weights
#' 
#' @export
#'
#' @examples
#' data_frame <- data.frame(x = simple_ts)
#' lags <- list(x = 0:1)
#' pbs_edm(data_frame, lags)
#' 
pbs_edm <- function(data_frame,
                    lags,
                    # from_user = seq_len(nrow(data_frame)), 
                    # into_user = seq_len(nrow(data_frame)), 
                    forecast_distance = 1L,
                    first_difference = FALSE,
                    centre_and_scale = FALSE,
                    show_calculations = FALSE) {
  
  # Check arguments
  assertthat::assert_that(
    is.data.frame(data_frame),
    is.list(lags),
    all(is.element(names(lags), names(data_frame))),
    is.count(forecast_distance),
    is.flag(first_difference),
    is.flag(centre_and_scale),
    is.flag(show_calculations)
  )
  
  # 
  nt_names <- names(lags)
  nt_matrix <- as.matrix(dplyr::select(tibble::as_tibble(data_frame), nt_names))
  
  xt_matrix <- nt_matrix
  
  if (first_difference) {
    xt_matrix <- rbind(apply(xt_matrix, 2, diff), NA_real_)
  }
  
  # Fix below here
  xt_mean <- apply(xt_matrix, 2, mean, na.rm = TRUE)
  xt_sd <- apply(xt_matrix, 2, sd, na.rm = TRUE)
  if (centre_and_scale) {
    xt_matrix <- apply(xt_matrix, 2, function(xt) (xt - xt_mean) / xt_sd)
  }
  
  # Calculate Xt from Nt
  nt_matrix <- as.matrix(data_frame[, data_names])

  
  
  nt <- dplyr::pull(data_frame, var = data_names[1])
  xt <- nt
  if (first_difference) {
    xt <- c(diff(xt), NA_real_)
  }
  if (centre_and_scale) {
    mean_xt <- mean(xt, na.rm = TRUE)
    sd_xt <- sd(xt, na.rm = TRUE)
    xt <- (xt - mean_xt) / sd_xt
  }
  xt_frame <- data.frame
  
  
  
  data_frame <- cbind(data.frame(Xt = xt), data_frame)
  
  
  trans_data <- 
    
    
    
  
  
  # Make lags matrix
  lag_sizes_vector <- unlist(lags, use.names = FALSE)
  col_names_vector <- rep(names(lags), lengths(lags))
  cols_list <- mapply(FUN = dplyr::pull,
                      var = col_names_vector,
                      MoreArgs = list(.data = data_frame),
                      SIMPLIFY = FALSE)
  lags_matrix <- mapply(FUN = dplyr::lag,
                        x = cols_list,
                        n = lag_sizes_vector)
  observations <- lags_matrix[, 1]
  
  # Make distance matrix
  lags_matrix[which(is.na(rowSums(lags_matrix))), ] <- NA # For distance
  distance_matrix <- as.matrix(dist(lags_matrix, diag = TRUE, upper = TRUE))
  diag(distance_matrix) <- NA
  
  # Exclude indices that project beyond the time series
  threshold <- ncol(distance_matrix) - forecast_distance
  distance_matrix[, which(seq_len(ncol(distance_matrix)) > threshold)] <- NA
  distance_matrix[which(seq_len(ncol(distance_matrix)) > threshold), ] <- NA
  
  # Exclude indices that contain the projection of the focal value
  # TODO: Check whether lags_unique is correct
  indices <- seq_len(nrow(distance_matrix))
  lags_unique <- unique(lag_sizes_vector)
  focal_indices <- rep(indices, each = length(lags_unique))
  exclude_indices <- focal_indices + forecast_distance + lags_unique
  within_range <- which(exclude_indices %in% indices)
  exclude_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                         y = exclude_indices[within_range]))
  distance_matrix[exclude_matrix] <- NA
  
  # Exclude indices that project to a row that contains NAs
  na_indices <- which(is.na(rowSums(lags_matrix)))
  exclude_indices <- na_indices - forecast_distance
  within_range <- which(exclude_indices %in% indices)
  focal_indices <- rep(indices, each = length(within_range))
  exclude_matrix <- as.matrix(data.frame(x = focal_indices,
                                         y = rep(exclude_indices[within_range],
                                                 length(focal_indices))))
  distance_matrix[exclude_matrix] <- NA
  
  # Exclude indices symmetrically around the forecast index
  if (symmetric_exclusion == TRUE) {
    symm_na_indices <- focal_indices + forecast_distance - lags_unique
    within_range <- which(symm_na_indices %in% indices)
    symm_na_matrix <- as.matrix(data.frame(x = focal_indices[within_range],
                                           y = symm_na_indices[within_range]))
    distance_matrix[symm_na_matrix] <- NA
  }

  # Make neighbour index matrix
  num_nbrs <- length(lags_unique) + 1
  seq_nbrs <- seq_len(num_nbrs)
  nbr_ind_matrix <- t(apply(distance_matrix, 1, order))[, seq_nbrs]
  nbr_ind_matrix[which(rowSums(!is.na(distance_matrix)) < num_nbrs), ] <- NA
  
  # Make neighbour value matrix
  nbr_val_matrix <- t(apply(nbr_ind_matrix, 1, function(x, y) y[x, 1],
                            y = lags_matrix)) # Works because NAs are NA_real_
  
  # Make neighbour distance matrix
  nbr_dst_matrix <- t(apply(distance_matrix, 1, sort, na.last = T))[, seq_nbrs]
  
  # Calculate weight matrix
  weight_matrix <- t(apply(nbr_dst_matrix, 1, function(x) exp(-x / x[1])))
  
  # Make projected neighbour index matrix
  pro_ind_matrix <- apply(nbr_ind_matrix, 2, dplyr::lag, n = forecast_distance)
  pro_ind_matrix <- pro_ind_matrix + forecast_distance
  
  # Make projected neighbour value matrix
  pro_val_matrix <- t(apply(pro_ind_matrix, 1, function(x, y) y[x, 1],
                            y = lags_matrix)) # Works because NAs are NA_real_
  
  # Make projected neighbour weight matrix
  pro_wt_matrix <- apply(weight_matrix, 2, dplyr::lag, n = forecast_distance)

  # Make forecasts
  forecasts <- rowSums(pro_val_matrix * pro_wt_matrix) / rowSums(pro_wt_matrix)
  forecasts <- as.vector(forecasts)
  
  # Compute statistics
  rho <- cor(observations, forecasts, use = "pairwise.complete.obs")
  rmse <- sqrt(mean((observations - forecasts)^2, na.rm = TRUE))
  
  # Results tibble
  results_only <- tibble::tibble(
    E = length(unlist(lags)),
    rho = rho,
    rmse = rmse)
  
  # Calculations list
  results_with_calculations <- list(
    results = results_only,
    lags = lags,
    Nt_observed,
    Nt_forecast = NULL,
    Xt_observed,
    Xt_forecast,
    lags_matrix = lags_matrix,
    nbr_index = nbr_ind_matrix,
    pro_index = pro_ind_matrix,
    nbr_distance = round(nbr_dst_matrix, 3),
    pro_distance = round(pro_dst_matrix, 3),
    nbr_weights = round(weight_matrix, 3),
    pro_weights = round(pro_wt_matrix, 3)
  )
  
  # Return value
  if (show_calculations) {
    results_with_calculations
  } else {
    results_only
  }
}
