#' Calculate Network Features
#'
#' This function calculates specific network features for the given igraph object.
#' The available features are: "clustering_coefficient", "global_efficiency",
#' "degrees", or "all".
#'
#' @param graph An igraph object representing the network. Can be created with `make_network()`.
#' @param feature A character string specifying the feature to compute.
#' Accepts one of: "clustering_coefficient", "global_efficiency",
#' "degrees", or "all". Defaults to "all".
#'
#' @return A list containing the requested feature(s) for the graph.
#' If "all" is selected, the list contains all three features.
#'
#' @examples
#' # Simulate a binarized calcium matrix
#' binarized_calcium_matrix <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#'
#' # Generate the network graph
#' graph <- make_network(binarized_calcium_matrix)
#'
#' # Calculate clustering coefficient
#' network_features(graph, feature = "clustering_coefficient")
#'
#' # Calculate global efficiency
#' network_features(graph, feature = "global_efficiency")
#'
#' # Calculate both
#' network_features(graph, feature = "all")
#'
#' @export
#' @importFrom igraph transitivity degree
#' @importFrom igraph global_efficiency
#'
network_features <- function(graph, feature = "all") {

  # Check if the graph is an igraph object
  if (!inherits(graph, "igraph")) {
    stop("The 'graph' argument must be an igraph object.")
  }

  # Initialize an empty list to store the features
  network_features <- list()

  # Calculate clustering coefficient if requested
  if (feature == "clustering_coefficient" || feature == "all") {
    clustcoeff <- igraph::transitivity(graph)
    network_features$clustering_coefficient <- clustcoeff
  }

  # Calculate global efficiency if requested
  if (feature == "global_efficiency" || feature == "all") {
    globaleff <- igraph::global_efficiency(graph)
    network_features$global_efficiency <- globaleff
  }

  # Calculate degrees if requested
  if (feature == "degrees" || feature == "all") {
    node_degrees <- igraph::degree(graph)
    network_features$degrees <- node_degrees
  }

  # If specific feature is requested, return it directly, else return all features
  if (feature != "all") {
    return(network_features[[feature]])
  } else {
    return(network_features)
  }
}
