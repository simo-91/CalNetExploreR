#' Subset and Analyze Connections Between Labeled and Unlabeled Cells
#'
#' This function calculates the number of connections between a subset of labeled cells and the rest of the unlabeled cells, based on a correlation matrix. It applies a correlation threshold to filter out weak connections and computes the total connections, possible connections, and their proportions for labeled-to-unlabeled and labeled-to-labeled cell pairs.
#'
#' @param correlation_matrix A square matrix representing the correlation between cells, where each entry contains the correlation value between two cells. The matrix should be of the same dimension as the number of cells in `coordinates`.
#' @param coordinates A dataframe containing the coordinates of the cells. It must contain a `Label` column, where labeled cells have non-zero values, and unlabeled cells have a value of 0.
#' @param correlation_threshold A numeric value to apply a threshold for correlations. Correlations below this value will be set to 0 (i.e., no connection), and correlations equal to or above this value will be set to 1 (i.e., a connection). Default is 0.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{total_connections_labeled_to_unlabeled}{The total number of connections from labeled cells to unlabeled cells.}
#'   \item{total_possible_connections_labeled_to_unlabeled}{The total possible number of connections from labeled to unlabeled cells.}
#'   \item{proportion_labeled_to_unlabeled}{The proportion of labeled-to-unlabeled connections relative to the total possible labeled-to-unlabeled connections.}
#'   \item{total_connections_labeled_to_labeled}{The total number of connections between labeled cells.}
#'   \item{total_possible_connections_labeled_to_labeled}{The total possible number of connections between labeled cells (excluding self-connections).}
#'   \item{proportion_labeled_to_labeled}{The proportion of labeled-to-labeled connections relative to the total possible labeled-to-labeled connections.}
#' }
#'
#' @examples
#' # Example correlation matrix (10x10) and corresponding coordinates
#' correlation_matrix <- matrix(runif(100, -1, 1), nrow = 10)
#' coordinates <- data.frame(
#'   X = runif(10),
#'   Y = runif(10),
#'   Label = c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0)
#' )
#'
#' # Apply subset_connections with a threshold of 0.5
#' result <- subset_connections(correlation_matrix, coordinates, correlation_threshold = 0.5)
#'
#' # View the result
#' print(result)
#' @export
subset_connections <- function(correlation_matrix, coordinates, correlation_threshold = 0) {
  # Check if the number of cells matches between the correlation matrix and coordinates
  if (nrow(correlation_matrix) != nrow(coordinates)) {
    stop("The number of cells in the correlation matrix and coordinates dataframe must be the same.")
  }

  # Extract labels and ensure they are in the correct order
  labels <- coordinates$Label
  if (is.null(labels)) {
    stop("The 'coordinates' dataframe must contain a 'Label' column.")
  }

  # Identify indices of labeled and unlabeled cells
  labeled_indices <- which(labels != 0)  # Assuming '0' indicates unlabeled cells
  unlabeled_indices <- which(labels == 0)

  # Apply the correlation threshold to create a binary adjacency matrix
  adj_matrix <- correlation_matrix
  adj_matrix[adj_matrix < correlation_threshold] <- 0
  adj_matrix[adj_matrix >= correlation_threshold] <- 1

  # Remove self-connections by setting the diagonal to zero
  diag(adj_matrix) <- 0

  # Calculate connections from labeled to unlabeled cells
  connections_labeled_to_unlabeled <- adj_matrix[labeled_indices, unlabeled_indices]
  total_connections_labeled_to_unlabeled <- sum(connections_labeled_to_unlabeled)

  # Calculate connections between labeled cells
  connections_labeled_to_labeled <- adj_matrix[labeled_indices, labeled_indices]
  total_connections_labeled_to_labeled <- sum(connections_labeled_to_labeled)

  # Calculate total possible connections
  num_labeled <- length(labeled_indices)
  num_unlabeled <- length(unlabeled_indices)
  total_possible_connections_labeled_to_unlabeled <- num_labeled * num_unlabeled
  total_possible_connections_labeled_to_labeled <- num_labeled * (num_labeled - 1)  # Excluding self-connections

  # Compute proportions
  proportion_labeled_to_unlabeled <- total_connections_labeled_to_unlabeled / total_possible_connections_labeled_to_unlabeled
  proportion_labeled_to_labeled <- total_connections_labeled_to_labeled / total_possible_connections_labeled_to_labeled

  # Return the results as a list
  return(list(
    total_connections_labeled_to_unlabeled = total_connections_labeled_to_unlabeled,
    total_possible_connections_labeled_to_unlabeled = total_possible_connections_labeled_to_unlabeled,
    proportion_labeled_to_unlabeled = proportion_labeled_to_unlabeled,
    total_connections_labeled_to_labeled = total_connections_labeled_to_labeled,
    total_possible_connections_labeled_to_labeled = total_possible_connections_labeled_to_labeled,
    proportion_labeled_to_labeled = proportion_labeled_to_labeled
  ))
}
