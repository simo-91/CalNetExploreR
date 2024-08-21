#' Generate Population Activity Plots
#'
#' Generates a raster plot and line plot after performing hierarchical clustering on the data to sort similar cells.
#'
#' @param binarized_calcium_matrix A binary matrix where each row represents a cell and each column represents a timepoint.
#' @param binarize A logical value indicating whether to binarize the calcium matrix. If TRUE, the matrix will be binarized using binarize(). Defaults to FALSE.
#' @param dendrogram A logical value indicating whether to include the dendrogram plot. Defaults to FALSE.
#' @return A combined plot showing the raster plot with hierarchical clustering and a line plot of population activity.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' plot <- population_activity.plt(calcium_matrix, binarize = TRUE, dendrogram = TRUE)
#' @export
#' @import ggplot2
#' @import reshape2
#' @import ggdendro
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
#' @importFrom cowplot align_plots plot_grid
#' @importFrom ggpubr theme_pubr
#' @importFrom grid viewport grid.newpage
population_activity.plt <- function(binarized_calcium_matrix, binarize = FALSE, dendrogram = FALSE) {
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
  hc <- stats::hclust(stats::dist(calcium_matrix_binary, method = "euclidean"), method = "ward.D2")
  dhc <- stats::as.dendrogram(hc)

  # Generate the dendrogram plot
  peaks.dendro <- ggdendro::ggdendrogram(dhc, rotate = TRUE, labels = FALSE) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  # Order the cells according to the hierarchical clustering
  peaks.order <- stats::order.dendrogram(dhc)
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
    ggplot2::ggtitle("Population activity over time")

  # Calculate the percentage of active cells over time
  total_cells <- nrow(calcium_matrix_binary)
  spksSUM <- colSums(calcium_matrix_binary) / total_cells * 100
  spksSUM_df <- data.frame(time = 0:(length(spksSUM) - 1), activity = spksSUM)

  # Plot the percentage of active cells over time
  spksSUM.plt <- ggplot2::ggplot(spksSUM_df, ggplot2::aes(time, activity)) +
    ggplot2::geom_line() +
    ggpubr::theme_pubr() +
    ggplot2::ylab("Active Cells %") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  # Combine the raster plot and the activity plot
  if (dendrogram) {
    grid::grid.newpage()
    print(raster.hc, vp = grid::viewport(x = 0.4, y = 0.6, width = 0.8, height = 0.8))
    print(peaks.dendro, vp = grid::viewport(x = 0.88, y = 0.58, width = 0.25, height = 0.86))
    print(spksSUM.plt, vp = grid::viewport(x = 0.4, y = 0.1, width = 0.8, height = 0.2))
  } else {
    plots <- cowplot::align_plots(raster.hc, spksSUM.plt, align = 'v', axis = 'l')
    combined_plot <- cowplot::plot_grid(plots[[1]], spksSUM.plt, ncol = 1, rel_heights = c(4.5, 1))
    return(combined_plot)
  }
}
