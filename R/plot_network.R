#' Plot a Network Graph from Calcium Imaging Data
#'
#' This function plots a network graph where nodes represent cells and edges represent connections between them.
#' The nodes can be color-coded based on the selected label, either by community membership or by the frequency of events per minute.
#'
#' @param graph An igraph object representing the network. Can be created with \code{make_network()}.
#' @param coordinates A data frame containing X and Y coordinates for each cell. Must include columns \code{"X"}, \code{"Y"}, and \code{"Cell"}.
#' @param label A character string indicating how to color-code the nodes. Options are \code{"communities"} or \code{"frequency"}. Defaults to \code{"communities"}.
#' @param cell_ID Either a character vector of labels for the nodes, or a special value to specify how to label the nodes. If set to \code{"none"}, the nodes will be labeled with their numbers. If set to \code{"communities"}, the nodes will be labeled with the community number they belong to. Defaults to \code{"none"}.
#' @param reverse_y_scale A logical value indicating whether to reverse the Y scale in the plot (useful for matching image coordinates). Defaults to \code{FALSE}.
#' @param frequency_values A numeric vector containing the frequency of events per minute for each cell. Required if \code{label = "frequency"}.
#' @param correlation_threshold A numeric value specifying the threshold for filtering edges by weight. Set to \code{"none"} to disable filtering. Defaults to \code{0.3}.
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
#'
#' # Plot the network graph with communities as labels inside nodes
#' plot <- plot_network(graph, coordinates = posXY, label = "communities", cell_ID = "communities")
#' print(plot)
#' @export
#' @importFrom igraph vcount delete_edges E cluster_leading_eigen
#' @importFrom grDevices colors heat.colors
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph scale_edge_color_viridis scale_edge_alpha_continuous
#' @importFrom ggplot2 aes scale_fill_manual scale_fill_gradientn scale_size_continuous theme margin unit element_text ggtitle
#' @importFrom RColorBrewer brewer.pal
plot_network <- function(graph, coordinates, label = "communities", cell_ID = "none", reverse_y_scale = FALSE, frequency_values = NULL, correlation_threshold = 0.3) {

  # Set layout for the graph using the provided coordinates
  layout <- as.matrix(coordinates[, c("X", "Y")])
  rownames(layout) <- coordinates$Cell

  # Apply the correlation threshold to filter edges, if specified
  if (correlation_threshold != "none") {
    graph <- igraph::delete_edges(graph, which(igraph::E(graph)$weight < correlation_threshold))
  }

  # Perform community detection if needed
  if (label == "communities" || cell_ID == "communities") {
    community_detection <- igraph::cluster_leading_eigen(graph, options = list(maxiter = 100000000))
    node_communities <- community_detection$membership
  }

  # Handle cell_ID assignment
  if (identical(cell_ID, "none")) {
    cell_ID <- as.character(1:igraph::vcount(graph))
  } else if (identical(cell_ID, "communities")) {
    cell_ID <- as.character(node_communities)
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

  # Determine how to color the nodes
  if (label == "communities") {
    node_colors <- factor(node_communities)

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
    ggraph::geom_node_text(ggplot2::aes(label = cell_ID), colour = "black", fontface = "plain", size = 3) +
    color_scale +
    ggplot2::scale_size_continuous(range = c(5, 12), guide = "none") +
    ggraph::theme_graph(background = "white", plot_margin = ggplot2::margin(5, 5, 5, 5)) +
    ggplot2::theme(legend.position = "right", legend.margin = ggplot2::margin(1, 1, 1, 1),
                   legend.key.size = ggplot2::unit(0.5, 'cm'),
                   text = ggplot2::element_text(family = "sans")) +
    ggplot2::ggtitle("Network Graph")

  # Reverse the Y scale if required
  if (reverse_y_scale) {
    graph.plt <- graph.plt + ggplot2::scale_y_reverse()
  }

  return(graph.plt)
}
