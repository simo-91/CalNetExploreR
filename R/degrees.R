#' Degree Distribution Analysis
#'
#' This function calculates the degree of each node in the network, generates a histogram of the degree distribution,
#' and returns either the plot or the mean degree based on the `plot` argument.
#'
#' @param graph An igraph object representing the network.
#' @param plot A logical value indicating whether to generate and return the degree histogram plot. Defaults to TRUE.
#' @return If `plot` is TRUE, a ggplot object representing the degree histogram. If `plot` is FALSE, a numeric value representing the mean degree.
#' @examples
#' binarized_calcium_matrix <- matrix(runif(100, 0, 1), nrow = 10, ncol = 10)
#' graph <- make_network(binarized_calcium_matrix)
#' result <- degree_analysis(graph, plot = TRUE)
#' @export
#' @importFrom igraph degree
#' @import ggplot2
#' @importFrom ggpubr gghistogram
degrees <- function(graph, plot = TRUE) {
  # Calculate degrees of nodes
  degrees <- igraph::degree(graph)

  if (plot) {
    # Create a histogram of degrees
    degree.hist <- ggpubr::gghistogram(data.frame(Degree = degrees), x = "Degree", y = "..count..", binwidth = 1) +
      ggplot2::ggtitle("Degree Distribution Histogram") +
      ggplot2::xlab("Degree") +
      ggplot2::ylab("Count") +
      ggplot2::theme_minimal()

    return(degree.hist)
  } else {
    # Calculate and return the mean degree
    degree.mean <- mean(degrees)
    return(degree.mean)
  }
}
