#' Binarize Calcium Imaging Data
#'
#' This function binarizes the timeseries calcium data for each cell using a specified cutoff function.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param cutoff_func A function to determine the threshold for binarizing the data. Default is twice the standard deviation of each cell.
#' @return A binary matrix where each cell's timeseries data is converted to 0 or 1 based on the cutoff function.
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' binary_data <- binarize(data)
#' @export
binarize <- function(calcium_matrix, cutoff_func = function(x) {
  th <- 2 * sd(x)
  x[x <= th] <- 0
  x[x > th] <- 1
  return(x)
}) {
  t(apply(calcium_matrix, 1, cutoff_func))
}
