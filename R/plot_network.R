#' Plot a Network Graph from Calcium Data
#'
#' This function plots a network graph where nodes represent cells and edges represent connections between them.
#' The nodes can be color-coded by community membership, and edges are color-coded by their weights.
#'
#' @param graph An igraph object representing the network. Can be created with `make_graph()`.
#' @param coordinates A data frame containing X and Y coordinates for each cell ID. Must include columns "X", "Y", and "Cell".
#' @param communities A logical value indicating whether to color-code nodes by community membership using `leading.eigenvector.community()`. Defaults to FALSE.
#' @param cell_ID A vector of cell IDs. Defaults to the row names of `calcium_matrix_binarized`.
#' @param reverse_y_scale A logical value indicating whether to reverse the Y scale in the plot (useful for matching image coordinates). Defaults to FALSE.
#' @return A ggplot object representing the network graph.
#' @examples
#' # Example usage:
#' graph <- make_graph(binarized_calcium_matrix)
#' posXY <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10)
#' plot <- plot_network(graph, coordinates = posXY, communities = TRUE)
#' print(plot)
#' @export
#' @import igraph
#' @import ggraph
#' @import ggplot2
#' @import viridis
#' @import RColorBrewer
plot_network <- function(graph, coordinates, communities = FALSE, cell_ID = rownames(calcium_matrix_binarized), reverse_y_scale = FALSE) {
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)

  # Ensure that the length of cell_ID matches the number of nodes in the graph
  if (length(cell_ID) != vcount(graph)) {
    stop("The length of cell_ID must match the number of nodes in the graph.")
  }

  # Set layout for the graph using the provided coordinates
  layout <- as.matrix(coordinates[, c("X", "Y")])
  rownames(layout) <- coordinates$Cell

  # Generate a large palette for community coloring
  colorz <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
  colorz.invert <- rev(colorz)
  communities.palette.big <- c(brewer.pal(8, 'Dark2'),
                               brewer.pal(8, 'Accent'),
                               brewer.pal(9,'Set1'),
                               brewer.pal(8,'Set2'),
                               brewer.pal(12,'Set3'),
                               brewer.pal(8, 'Accent'),
                               brewer.pal(9,'Pastel1'),
                               colorz.invert, brewer.pal(8, 'Accent'),
                               colorz)

  # Calculate community membership if required
  if (communities) {
    community_detection <- leading.eigenvector.community(graph, options = list(maxiter = 100000000))
    membership <- factor(community_detection$membership)

    # Calculate community sizes
    community_sizes <- table(membership)

    # Filter out smaller communities (e.g., less than 5 members)
    large_communities <- names(community_sizes[community_sizes >= 2])
    membership <- factor(membership, levels = large_communities)
  } else {
    membership <- NULL
  }

  # Create the network plot using ggraph
  graph.plt <- ggraph(graph, layout = layout) +
    geom_edge_link(aes(colour = weight, alpha = weight)) +
    scale_edge_alpha_continuous(range = c(0.1, 1), guide = "none") +
    scale_edge_color_viridis(name = "F. Corr", alpha = 1, begin = 0.3, end = 1,
                             discrete = FALSE, option = "inferno", direction = 1,
                             guide = guide_colourbar(available_aes = "edge_colour")) +
    geom_node_point(aes(fill = membership, size = degree(graph)), shape = 21) +
    geom_node_text(aes(label = cell_ID), colour = "black", fontface = 1, size = 3) +
    scale_fill_manual(values = communities.palette.big, drop = FALSE) +
    scale_size_continuous(range = c(5, 12), guide = "none") +
    theme_graph(background = "white", plot_margin = margin(5, 5, 5, 5)) +
    theme(legend.position = "right", legend.margin = margin(1, 1, 1, 1),
          legend.key.size = unit(0.5, 'cm')) +
    ggtitle("Network Graph")

  # Reverse the Y scale if required
  if (reverse_y_scale) {
    graph.plt <- graph.plt + scale_y_reverse()
  }

  # Add a guide for the larger communities
  if (communities) {
    graph.plt <- graph.plt + guides(fill = guide_legend(title = "Communities",
                                                        override.aes = list(size = 5),
                                                        ncol = 1))
  }

  return(graph.plt)
}
