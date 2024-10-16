#' Calculate Normalized Labeled-to-NonLabeled Connections
#'
#' This function calculates the number of Labeled-to-NonLabeled connections (where labeled nodes have any non-zero value in the Label column)
#' and normalizes the result by the mean degree of the network.
#'
#' @param network An igraph object representing the network.
#' @param coordinates A data frame containing the coordinates and labeling information of the nodes. Must include a "Label" and "Cell" column.
#' @return A numeric value representing the normalized Labeled-to-NonLabeled connections.
#' @export
#' @importFrom igraph degree vcount

labeled_to_nonlabeled_connections_normalized <- function(network, coordinates) {
  # Ensure that the coordinates dataframe has a "Label" and "Cell" column
  if (!all(c("Label", "Cell") %in% colnames(coordinates))) {
    stop("'Label' and 'Cell' columns are missing from the coordinates dataframe.")
  }

  # Ensure that the number of vertices in the graph matches the number of rows in coordinates
  if (igraph::vcount(network) != nrow(coordinates)) {
    stop("Mismatch between the number of vertices in the graph and the number of rows in coordinates.")
  }

  # Calculate the mean degree of the network
  mean_degree <- mean(igraph::degree(network))

  # Get the total number of labeled-to-nonlabeled connections
  total_labeled_to_nonlabeled <- labeled_to_nonlabeled_connections(network, coordinates)

  # Normalize by the mean degree
  normalized_connections <- total_labeled_to_nonlabeled / mean_degree

  # Return the normalized connections
  return(normalized_connections)
}
