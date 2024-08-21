#' Create a Network object from a binarized calcium matrix
#'
#' This function produces a network based on the maximum cross-correlation between cells' calcium activity.
#' The user can specify the lag for the cross-correlation function and the correlation threshold for filtering edges.
#'
#' @param binarized_calcium_matrix A binarized matrix where each row represents a cell and each column represents a timepoint. This matrix can be generated using the `binarize()` function.
#' @param lag.max The maximum lag to use in the cross-correlation function (CCF). Defaults to 1.
#' @param correlation_threshold The threshold value for filtering edges in the network (Pearson's coefficients go from -1 to +1). Set to "none" to disable filtering. Defaults to "none".
#' @return An `igraph` object representing the network of correlated cells.
#' @examples
#' binarized_calcium_matrix <- matrix(sample(c(0, 1), 1000, replace = TRUE), nrow = 10)
#' network <- make_network(binarized_calcium_matrix)
#' network_no_filter <- make_network(binarized_calcium_matrix, correlation_threshold = "none")
#' @export
#' @importFrom igraph graph.adjacency delete.edges E
#' @importFrom utils txtProgressBar setTxtProgressBar
make_network <- function(binarized_calcium_matrix, lag.max = 1, correlation_threshold = "none") {

  # Transpose the matrix for cross-correlation calculation
  T.allcellst <- t(binarized_calcium_matrix)

  # Define a function to find max CCF between two signals within a specified lag range
  max_CCF <- function(a, b, lag.max) {
    d <- stats::ccf(a, b, plot = FALSE, lag.max = lag.max)
    cor <- d$acf[,,1]
    return(max(cor))
  }

  # Initialize the result matrix
  n_cols <- ncol(T.allcellst)
  cmat.allcellst <- matrix(0, n_cols, n_cols)

  # Sequential computation
  pb <- utils::txtProgressBar(min = 0, max = n_cols, style = 3)
  for (i in 1:(n_cols - 1)) {
    for (j in (i + 1):n_cols) {
      cmat.allcellst[i, j] <- max_CCF(T.allcellst[, i], T.allcellst[, j], lag.max)
    }
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  # Mirror the upper triangular part to the lower triangular part
  cmat.allcellst[lower.tri(cmat.allcellst)] <- t(cmat.allcellst)[lower.tri(cmat.allcellst)]
  diag(cmat.allcellst) <- 1  # Set diagonal to 1

  # Replace NaNs with 0 (to handle cases where the time-series might be constant)
  cmat.allcellst[is.na(cmat.allcellst)] <- 0

  # Create a network from the correlation matrix
  network <- igraph::graph.adjacency(as.matrix(cmat.allcellst), mode = "undirected", weighted = TRUE, diag = FALSE)

  # Apply the correlation threshold to filter edges, if specified
  if (correlation_threshold != "none") {
    network <- igraph::delete.edges(network, which(igraph::E(network)$weight < correlation_threshold))
  }

  return(network)
}
