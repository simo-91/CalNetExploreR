#' Calculate Percentage of Active Cells Over Time
#'
#' This function calculates the percentage of active cells over time from a binarized calcium matrix.
#'
#' @param binarized_calcium_matrix A binarized matrix where each row represents a cell and each column represents a timepoint. Can be created using binarize()
#' @return A data frame with the sum of active cells at each time point, the corresponding time, and the percentage of active cells.
#' @examples
#' binarized_data <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#' active_percentage <- active_cells_percentage(binarized_data)
#' @export
active_cells_percentage <- function(binarized_calcium_matrix) {
  active_cells <- colSums(binarized_calcium_matrix)
  active_cells_df <- as.data.frame(active_cells)
  active_cells_df$Time <- 0:(nrow(active_cells_df) - 1)
  active_cells_df$Perc <- active_cells_df$active_cells / nrow(binarized_calcium_matrix) * 100
  return(active_cells_df)
}
