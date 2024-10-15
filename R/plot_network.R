#' Plot a Network Graph from Calcium Data
#'
#' This function plots a network graph where nodes represent cells and edges represent connections between them.
#' The nodes can be color-coded based on the selected label, either by community membership or by the frequency of events per minute.
#'
#' @param graph An igraph object representing the network. Can be created with `make_network()`.
#' @param coordinates A data frame containing X and Y coordinates for each cell ID. Must include columns "X", "Y", and "Cell".
#' @param label A character string indicating what to label the cells with. Options are "communities" or "frequency". Defaults to "communities".
#' @param cell_ID A dataframe of cell IDs (should contain X and Y columns). If set to "none", the nodes will be labeled with their numbers. Defaults to "none".
#' @param reverse_y_scale A logical value indicating whether to reverse the Y scale in the plot (useful for matching image coordinates). Defaults to FALSE.
#' @param frequency_values A numeric vector containing the frequency of events per minute for each cell. Required if `label = "frequency"`.
#' @param correlation_threshold A numeric value specifying the threshold for filtering edges by weight. Set to "none" to disable filtering. Defaults to 0.3.
#' @return A ggplot object representing the network graph.
#' @examples
#' # Simulate a binarized calcium matrix
#' binarized_calcium_matrix <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#'
#' # Generate the network graph
#' graph <- make_network(binarized_calcium_matrix)
#'
#' # Simulate XY coordinates for the cells
#' posXY <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10)
#'
#' # Simulate frequency values for the cells
#' frequency_values <- runif(10, 0, 5)
#'
#' # Plot the network graph with frequency as the label
#' plot <- plot_network(graph, coordinates = posXY, label = "frequency", frequency_values = frequency_values)
#' print(plot)
#' @export
#' @importFrom igraph vcount delete_edges E
#' @importFrom grDevices colors
#' @importFrom grDevices heat.colors
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph scale_edge_color_viridis scale_edge_alpha_continuous
#' @importFrom ggplot2 aes scale_fill_manual scale_fill_gradientn scale_size_continuous theme margin unit element_text
#' @importFrom ggplot2 ggtitle
#' @importFrom RColorBrewer brewer.pal
plot_network <- function(graph, coordinates, label = "communities", cell_ID = "none", reverse_y_scale = FALSE, frequency_values = NULL, correlation_threshold = 0.3) {

  # Set layout for the graph using the provided coordinates
  layout <- as.matrix(coordinates[, c("X", "Y")])
  rownames(layout) <- coordinates$Cell

  # Apply the correlation threshold to filter edges, if specified
  if (correlation_threshold != "none") {
    graph <- igraph::delete_edges(graph, which(igraph::E(graph)$weight < correlation_threshold))
  }

  # Handle cell_ID assignment
  if (identical(cell_ID, "none")) {
    cell_ID <- as.character(1:igraph::vcount(graph))
  } else if (length(cell_ID) != igraph::vcount(graph)) {
    stop("The length of cell_ID must match the number of nodes in the graph.")
  }

  # Generate a large palette for community coloring (only used if label = "communities")
  colorz <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
  colorz.invert <- rev(colorz)
  communities.palette.big <- c(RColorBrewer::brewer.pal(8, 'Dark2'),
                               RColorBrewer::brewer.pal(9,'Set1'),
                               RColorBrewer::brewer.pal(8,'Set2'),
                               RColorBrewer::brewer.pal(12,'Set3'),
                               RColorBrewer::brewer.pal(8, 'Accent'),
                               colorz,
                               RColorBrewer::brewer.pal(9,'Pastel1'),
                               colorz.invert, colorz)

  # Determine how to label the nodes
  if (label == "communities") {
    community_detection <- igraph::cluster_leading_eigen(graph, options = list(maxiter = 100000000))
    node_colors <- factor(community_detection$membership)

    # Calculate community sizes
    community_sizes <- table(node_colors)
    community_sizes_numeric <- as.numeric(community_sizes[as.character(node_colors)])

    # Filter out smaller communities (e.g., less than 2 members)
    large_communities <- names(community_sizes[community_sizes >= 2])
    node_colors <- factor(node_colors, levels = large_communities)

    # Assign node sizes based on community size
    node_sizes <- community_sizes_numeric

    color_scale <- ggplot2::scale_fill_manual(values = communities.palette.big, drop = FALSE, guide = ggplot2::guide_legend(title = "Communities", ncol = 1))

    # Use "inferno" color scale for edges if communities are labeled
    edge_color_scale <- ggraph::scale_edge_color_viridis(name = "F. Corr", alpha = 1, begin = 0.3, end = 1,
                                                         discrete = FALSE, option = "inferno", direction = -1,
                                                         guide = ggraph::guide_edge_colourbar(available_aes = "edge_colour"))
  } else if (label == "frequency") {
    if (is.null(frequency_values)) {
      stop("You must provide 'frequency_values' when label is set to 'frequency'.")
    }
    if (length(frequency_values) != igraph::vcount(graph)) {
      stop("The length of 'frequency_values' must match the number of nodes in the graph.")
    }

    node_colors <- frequency_values
    node_sizes <- frequency_values
    color_scale <- ggplot2::scale_fill_gradientn(name = "Frequency (events/min)", colors = rev(grDevices::heat.colors(5)))

    # Use "mako" color scale for edges if frequency is labeled
    edge_color_scale <- ggraph::scale_edge_color_viridis(name = "F. Corr", alpha = 1, begin = 0.3, end = 1,
                                                         discrete = FALSE, option = "mako", direction = -1,
                                                         guide = ggraph::guide_edge_colourbar(available_aes = "edge_colour"))
  } else {
    stop("Invalid 'label' option. Choose either 'communities' or 'frequency'.")
  }

  # Create the network plot using ggraph
  graph.plt <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(ggplot2::aes(colour = weight, alpha = weight)) +
    ggraph::scale_edge_alpha_continuous(range = c(0.1, 1), guide = "none") +
    edge_color_scale +
    ggraph::geom_node_point(ggplot2::aes(fill = node_colors, size = node_sizes), shape = 21) +
    ggraph::geom_node_text(ggplot2::aes(label = cell_ID), colour = "black", fontface = "plain", size = 3) +  # Use a common font
    color_scale +
    ggplot2::scale_size_continuous(range = c(5, 12), guide = "none") +
    ggraph::theme_graph(background = "white", plot_margin = ggplot2::margin(5, 5, 5, 5)) +
    ggplot2::theme(legend.position = "right", legend.margin = ggplot2::margin(1, 1, 1, 1),
                   legend.key.size = ggplot2::unit(0.5, 'cm'),
                   text = ggplot2::element_text(family = "sans")) +  # Set default font family
    ggplot2::ggtitle("Network Graph")

  # Reverse the Y scale if required
  if (reverse_y_scale) {
    graph.plt <- graph.plt + ggplot2::scale_y_reverse()
  }

  return(graph.plt)
}
