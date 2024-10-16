#' Calculate Labeled-to-NonLabeled Connections
#'
#' This function calculates the total number of edges that connect labeled nodes (any non-zero value in the Label column) to non-labeled nodes (Label == 0).
#'
#' @param network An igraph object representing the network.
#' @param coordinates A data frame containing the coordinates and labeling information of the nodes. Must include a "Label" and "Cell" column.
#' @return A numeric value representing the total number of Labeled-to-NonLabeled connections.
#' @export
#' @importFrom igraph neighbors vcount

labeled_to_nonlabeled_connections <- function(network, coordinates) {
  # Ensure that the coordinates dataframe has a "Label" and "Cell" column
  if (!all(c("Label", "Cell") %in% colnames(coordinates))) {
    stop("'Label' and 'Cell' columns are missing from the coordinates dataframe.")
  }

  # Ensure that the number of vertices in the graph matches the number of rows in coordinates
  if (igraph::vcount(network) != nrow(coordinates)) {
    stop("Mismatch between the number of vertices in the graph and the number of rows in coordinates.")
  }

  # Identify the indices of labeled and non-labeled nodes
  labeled_indices <- which(coordinates$Label != 0)  # Any non-zero value is considered labeled
  non_labeled_indices <- which(coordinates$Label == 0)

  # Initialize a counter for connections
  labeled_to_non_labeled_connections <- 0

  # Iterate over each labeled node and count its connections to non-labeled nodes
  for (labeled_index in labeled_indices) {
    # Get the corresponding vertex ID from the Cell column (adjust for 0-based indexing in coordinates and 1-based in igraph)
    labeled_vertex <- coordinates$Cell[labeled_index] + 1

    # Ensure the labeled vertex exists in the graph
    if (labeled_vertex > igraph::vcount(network)) {
      stop("Invalid labeled_vertex index for graph.")
    }

    # Get the neighbors of the labeled node
    neighbors <- igraph::neighbors(network, labeled_vertex, mode = "all")

    # Convert neighbors to their corresponding Cell IDs (adjust for 0-based indexing in coordinates)
    neighbor_indices <- as.numeric(neighbors) - 1

    # Count the connections to non-labeled nodes
    labeled_to_non_labeled_connections <- labeled_to_non_labeled_connections + sum(neighbor_indices %in% coordinates$Cell[non_labeled_indices])
  }

  # Return the total number of Labeled-to-NonLabeled connections
  return(labeled_to_non_labeled_connections)
}
