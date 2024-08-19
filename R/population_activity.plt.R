#' Generate Population Activity Plots
#'
#' Generates a raster plot and line plot after performing hierarchical clustering on the data to sort similar cells.
#'
#' @param binarized_calcium_matrix A binary matrix where each row represents a cell and each column represents a timepoint.
#' @param binarize A logical value indicating whether to binarize the calcium matrix. If TRUE, the matrix will be binarized using binarize(). Defaults to FALSE.
#' @return A combined plot showing the raster plot with hierarchical clustering and a line plot of population activity.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' plot <- population_activity.plt(calcium_matrix, binarize = TRUE)
#' @export
#' @import ggplot2
#' @import reshape2
#' @import ggdendro
#' @import cowplot
#' @import ggpubr
population_activity.plt <- function(binarized_calcium_matrix, binarize = FALSE) {
  # If binarize is TRUE, binarize the calcium matrix
  if (binarize) {
    calcium_matrix_binary <- binarize(binarized_calcium_matrix)
  } else {
    calcium_matrix_binary <- binarized_calcium_matrix
  }

  # Transform the calcium matrix for plotting
  dfpeaks <- as.data.frame(t(calcium_matrix_binary))
  dfpeaks$time <- 0:(nrow(dfpeaks) - 1)
  meltPeaks <- reshape2::melt(dfpeaks, id = "time")
  colnames(meltPeaks) <- c('time', 'cell', 'Ca2+')

  # Perform hierarchical clustering
  hc <- hclust(dist(calcium_matrix_binary, method = "euclidean"), method = "ward.D2")
  dhc <- as.dendrogram(hc)

  # Generate the dendrogram plot
  peaks.dendro <- ggdendro::ggdendrogram(dhc, rotate = TRUE, labels = FALSE) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  # Order the cells according to the hierarchical clustering
  peaks.order <- order.dendrogram(dhc)
  peaks.rows <- rownames(calcium_matrix_binary)
  peaks.rows <- as.data.frame(peaks.rows)
  meltPeaks$cell <- factor(x = meltPeaks$cell,
                           levels = peaks.rows$peaks.rows[peaks.order],
                           ordered = TRUE)

  # Generate the raster plot
  raster.hc <- ggplot2::ggplot(meltPeaks, ggplot2::aes(time, cell)) +
    ggplot2::geom_raster(ggplot2::aes(fill = `Ca2+`)) +
    ggplot2::scale_fill_gradientn(colours = c("white", "black")) +
    ggpubr::theme_pubr() +
    ggplot2::ylab("cell ID") +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(colour = "red", hjust = .5)) +
    ggplot2::ggtitle(paste0("Population activity over time"))

  # Calculate the sum of active cells over time
  spksSUM <- colSums(calcium_matrix_binary > 0)
  spksSUM_df <- data.frame(time = 0:(length(spksSUM) - 1), activity = spksSUM)
  spksSUM.plt <- ggplot2::ggplot(spksSUM_df, ggplot2::aes(time, activity)) +
    ggplot2::geom_line() +
    ggpubr::theme_pubr() +
    ggplot2::ylab("Active cells") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  # Combine the raster plot and the activity plot
  plots <- cowplot::align_plots(raster.hc, spksSUM.plt, align = 'v', axis = 'l')
  combined_plot <- cowplot::plot_grid(plots[[1]], spksSUM.plt, ncol = 1, rel_heights = c(4.5, 1))

  return(combined_plot)
}
