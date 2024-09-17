#' Get Subset of Cells Based on Label
#'
#' This function checks if the provided coordinates dataframe contains a "Label" column.
#' The Label column should contain either a combination of 0 and 1 (for FALSE and TRUE) or logical values (FALSE and TRUE).
#' Cells marked as TRUE or 1 are considered part of the subset.
#'
#' @param coordinates A dataframe containing X, Y, Cell ID and a Label columns.
#' @return A dataframe containing only the rows where Label is 1 or TRUE. If the Label column does not exist or is not properly formatted, an error is raised.
#' @examples
#' # Example coordinates with a Label column
#' coordinates <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10, Label = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0))
#' subset_cells <- get_subset(coordinates)
#' print(subset_cells)
#' @export
get_subset <- function(coordinates) {
  # Check if the "Label" column exists
  if (!"Label" %in% colnames(coordinates)) {
    stop("The 'Label' column is missing from the coordinates dataframe.")
  }

  # Check if Label contains appropriate values (0/1 or TRUE/FALSE)
  if (!all(coordinates$Label %in% c(0, 1, TRUE, FALSE))) {
    stop("The 'Label' column must contain only 0, 1, TRUE, or FALSE values.")
  }

  # Filter the rows where Label is 1 or TRUE
  subset_df <- coordinates[coordinates$Label %in% c(1, TRUE), ]

  # Return the subset dataframe
  return(subset_df)
}
