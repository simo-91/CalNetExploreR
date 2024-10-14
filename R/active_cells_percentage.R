#' Calculate the Percentage of Active Cells Over Time
#'
#' This function calculates the percentage of active cells over time in a binarized calcium matrix.
#' It can optionally plot the percentage of active cells over time.
#'
#' @param calcium_matrix_binarized A binarized matrix where each row represents a cell and each column represents a timepoint. This matrix can be generated using the `binarize()` function.
#' @param binarize A logical value indicating whether to binarize the calcium matrix. If TRUE, the function will apply the `binarize()` function to the calcium_matrix_binarized before calculation. Defaults to FALSE.
#' @param plot A logical value indicating whether to generate a plot of the percentage of active cells over time. Defaults to FALSE.
#' @return If `plot` is FALSE, returns a data frame containing the time points and the percentage of active cells. If `plot` is TRUE, returns a ggplot object of the percentage of active cells over time.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' result <- active_cells_percentage(calcium_matrix, binarize = TRUE)
#' plot <- active_cells_percentage(calcium_matrix, binarize = TRUE, plot = TRUE)
#' @export
#' @import ggplot2
#' @import ggpubr
active_cells_percentage <- function(calcium_matrix_binarized, binarize = FALSE, plot = FALSE) {
  # Optionally binarize the matrix
  if (binarize) {
    calcium_matrix_binarized <- binarize(calcium_matrix_binarized)
  }

  # Calculate the percentage of active cells over time
  total_cells <- nrow(calcium_matrix_binarized)
  spksSUM <- colSums(calcium_matrix_binarized) / total_cells * 100
  spksSUM_df <- data.frame(time = 0:(length(spksSUM) - 1), Coactive_Cells = spksSUM)

  # If plot is TRUE, generate and return the plot
  if (plot) {
    spksSUM.plt <- ggplot2::ggplot(spksSUM_df, ggplot2::aes(time, Coactive_Cells)) +
      ggplot2::geom_line() +
      ggpubr::theme_pubr() +
      ggplot2::ylab("Active Cells %") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
    return(spksSUM.plt)
  } else {
    return(spksSUM_df)
  }
}
