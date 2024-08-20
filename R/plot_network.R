#' Plot a Network Graph from Calcium Data
#'
#' This function plots a network graph where nodes represent cells and edges represent connections between them.
#' The nodes can be color-coded based on the selected label, either by community membership or by the frequency of events per minute.
#'
#' @param graph An igraph object representing the network. Can be created with `make_graph()`.
#' @param coordinates A data frame containing X and Y coordinates for each cell ID. Must include columns "X", "Y", and "Cell".
#' @param label A character string indicating what to label the cells with. Options are "communities" or "frequency". Defaults to "communities".
#' @param cell_ID A vector of cell IDs. Defaults to the row names of `calcium_matrix_binarized`.
#' @param reverse_y_scale A logical value indicating whether to reverse the Y scale in the plot (useful for matching image coordinates). Defaults to FALSE.
#' @param frequency_values A numeric vector containing the frequency of events per minute for each cell. Required if `label = "frequency"`.
#' @return A ggplot object representing the network graph.
#' @examples
#' # Example usage:
#' graph <- make_graph(binarized_calcium_matrix)
#' posXY <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10)
#' frequency_values <- runif(10, 0, 5) # Simulated frequency data
#' plot <- plot_network(graph, coordinates = posXY, label = "frequency", frequency_values = frequency_values)
#' print(plot)
#' @export
#' @import igraph
#' @import ggraph
#' @import ggplot2
#' @import viridis
#' @import RColorBrewer
plot_network <- function(graph, coordinates, label = "communities", cell_ID = rownames(calcium_matrix_binarized), reverse_y_scale = FALSE, frequency_values = NULL) {
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

  # Generate a large palette for community coloring (only used if label = "communities")
  colorz <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
  colorz.invert <- rev(colorz)
  communities.palette.big <- c(brewer.pal(8, 'Dark2'),
                               brewer.pal(9,'Set1'),
                               brewer.pal(8,'Set2'),
                               brewer.pal(12,'Set3'),
                               brewer.pal(8, 'Accent'),
                               colorz,
                               brewer.pal(9,'Pastel1'),
                               colorz.invert, colorz)

  # Determine how to label the nodes
  if (label == "communities") {
    community_detection <- leading.eigenvector.community(graph, options = list(maxiter = 100000000))
    node_colors <- factor(community_detection$membership)

    # Calculate community sizes
    community_sizes <- table(node_colors)
    community_sizes_numeric <- as.numeric(community_sizes[as.character(node_colors)])

    # Filter out smaller communities (e.g., less than 2 members)
    large_communities <- names(community_sizes[community_sizes >= 2])
    node_colors <- factor(node_colors, levels = large_communities)

    # Assign node sizes based on community size
    node_sizes <- community_sizes_numeric

    color_scale <- scale_fill_manual(values = communities.palette.big, drop = FALSE, guide = guide_legend(title = "Communities", ncol = 1))

    # Use "inferno" color scale for edges if communities are labeled
    edge_color_scale <- scale_edge_color_viridis(name = "F. Corr", alpha = 1, begin = 0.3, end = 1,
                                                 discrete = FALSE, option = "inferno", direction = -1,
                                                 guide = guide_colourbar(available_aes = "edge_colour"))
  } else if (label == "frequency") {
    if (is.null(frequency_values)) {
      stop("You must provide 'frequency_values' when label is set to 'frequency'.")
    }
    if (length(frequency_values) != vcount(graph)) {
      stop("The length of 'frequency_values' must match the number of nodes in the graph.")
    }

    node_colors <- frequency_values
    node_sizes <- frequency_values
    color_scale <- scale_fill_gradientn(name = "Frequency (events/min)", colors = rev(heat.colors(5)))

    # Use "mako" color scale for edges if frequency is labeled
    edge_color_scale <- scale_edge_color_viridis(name = "F. Corr", alpha = 1, begin = 0.3, end = 1,
                                                 discrete = FALSE, option = "mako", direction = -1,
                                                 guide = guide_colourbar(available_aes = "edge_colour"))
  } else {
    stop("Invalid 'label' option. Choose either 'communities' or 'frequency'.")
  }

  # Create the network plot using ggraph
  graph.plt <- ggraph(graph, layout = layout) +
    geom_edge_link(aes(colour = weight, alpha = weight)) +
    scale_edge_alpha_continuous(range = c(0.1, 1), guide = "none") +
    edge_color_scale +
    geom_node_point(aes(fill = node_colors, size = node_sizes), shape = 21) +
    geom_node_text(aes(label = cell_ID), colour = "black", fontface = 1, size = 3) +
    color_scale +
    scale_size_continuous(range = c(5, 12), guide = "none") +
    theme_graph(background = "white", plot_margin = margin(5, 5, 5, 5)) +
    theme(legend.position = "right", legend.margin = margin(1, 1, 1, 1),
          legend.key.size = unit(0.5, 'cm')) +
    ggtitle("Network Graph")

  # Reverse the Y scale if required
  if (reverse_y_scale) {
    graph.plt <- graph.plt + scale_y_reverse()
  }

  return(graph.plt)
}
