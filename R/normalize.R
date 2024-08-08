#' Normalize Calcium Imaging Data
#'
#' This function normalizes the timeseries calcium data for each cell.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @return A normalized matrix where each cell's timeseries data is scaled to [0, 1].
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' normalized_data <- normalize(data)
#' @export
normalize <- function(calcium_matrix) {
  t(apply(calcium_matrix, 1, function(x) (x - min(x)) / (max(x) - min(x))))
}
