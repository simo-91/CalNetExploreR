#' Calculate and Plot Percentage of Coactive Cells Over Time
#'
#' This function calculates the percentage of coactive cells over time from a binarized calcium matrix.
#' It returns a data frame with the time points and corresponding percentage of active cells. Optionally, it can also generate a plot.
#'
#' @param binarized_calcium_matrix A binary matrix where each row represents a cell and each column represents a timepoint.
#' @param binarize A logical value indicating whether to binarize the calcium matrix. If TRUE, the matrix will be binarized using binarize(). Defaults to FALSE.
#' @param plot A logical value indicating whether to generate a plot of the percentage of coactive cells over time. Defaults to FALSE.
#' @return A data frame showing the percentage of coactive cells at each timepoint.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' coactive_cells.df <- coactive_cells(calcium_matrix, binarize = TRUE, plot = TRUE)
#' @export
#' @import ggplot2
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
coactive_cells <- function(binarized_calcium_matrix, binarize = FALSE, plot = FALSE) {
  # If binarize is TRUE, binarize the calcium matrix
  if (binarize) {
    calcium_matrix_binary <- binarize(binarized_calcium_matrix)
  } else {
    calcium_matrix_binary <- binarized_calcium_matrix
  }

  # Calculate the percentage of coactive cells over time
  total_cells <- nrow(calcium_matrix_binary)
  spksSUM <- colSums(calcium_matrix_binary) / total_cells * 100
  coactive_cells.df <- data.frame(time = 0:(length(spksSUM) - 1), activity = spksSUM)

  # Plot if requested
  if (plot) {
    p <- ggplot2::ggplot(coactive_cells.df, ggplot2::aes(time, activity)) +
      ggplot2::geom_line() +
      ggpubr::theme_pubr() +
      ggplot2::ylab("Active Cells %") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
      ggplot2::ggtitle("Percentage of Coactive Cells Over Time")

    print(p)
  }

  return(coactive_cells.df)
}
