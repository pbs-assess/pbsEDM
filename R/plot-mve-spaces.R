#' Plot the axes used in the top multiview embedding reconstructions
#'
#' Takes the results list from [multiview_embedding()] and indicates the
#' variables and lags used in the top reconstruction subsets.
#'
#' @param res_mve list object resulting from a call to [multiview_embedding()]
#' @param col character; color to use for tiles where the lag is `TRUE`. Default "blue".
#' @param rho_on_top logical; if `TRUE` (default), rho values are along the top (x-axis).
#'   If `FALSE`, rho values are on the y-axis.
#' @param ... passed onto [create_named_lags()], namely `order` to order the
#' variables consistently.
#' @examples
#' \dontrun{
#'   plot_mve_spaces(age11_res)
#'   plot_mve_spaces(age11_res, rho_on_top = FALSE)
#' }
#'
#' @author Andrew M. Edwards and GitHub Copilot
plot_mve_spaces <- function(res_mve,
                            col = "blue",
                            rho_on_top = TRUE,
                            ...){
  tib <- make_mve_spaces_tibble(res_mve,
                         ...)

  # Reshape tibble to long format for ggplot
  tib_long <- tib %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -c(rho, row_id),
                        names_to = "lag_name",
                        values_to = "present")

  # Create mapping of row_id to rounded rho values for labels
  rho_labels <- tib %>%
    dplyr::mutate(row_id = dplyr::row_number(),
                  rho_label = sprintf("%.2f", rho)) %>%
    dplyr::select(row_id, rho_label) %>%
    tibble::deframe()

  if (rho_on_top) {
    # Rho values on top (x-axis), lag variables on y-axis
    p <- ggplot2::ggplot(tib_long,
                         ggplot2::aes(x = factor(row_id, levels = unique(row_id)),
                                     y = factor(lag_name, levels = rev(unique(lag_name))))) +
      ggplot2::geom_tile(ggplot2::aes(fill = present),
                         color = "white",
                         linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c("TRUE" = col, "FALSE" = "white"),
                                guide = "none") +
      ggplot2::scale_x_discrete(position = "top",
                               labels = function(x) rho_labels[as.numeric(x)]) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(x = "Rho",
                    y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0.5))
  } else {
    # Rho values on y-axis, lag variables on top (x-axis) - original approach
    p <- ggplot2::ggplot(tib_long,
                         ggplot2::aes(x = lag_name,
                                     y = factor(row_id, levels = rev(unique(row_id))))) +
      ggplot2::geom_tile(ggplot2::aes(fill = present),
                         color = "white",
                         linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c("TRUE" = col, "FALSE" = "white"),
                                guide = "none") +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::scale_y_discrete(labels = function(x) rho_labels[as.numeric(x)]) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(x = "",
                    y = "Rho") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0.5))
  }

  return(p)
}
