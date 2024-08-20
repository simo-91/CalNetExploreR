#' Create a Network object from a binarized calcium matrix
#'
#' This function produces a network based on the maximum cross-correlation between cells' calcium activity.
#' It supports both parallel and sequential processing. The user can specify the lag for the cross-correlation
#' function and the correlation threshold for filtering edges.
#'
#' @param binarized_calcium_matrix A binarized matrix where each row represents a cell and each column represents a timepoint. This matrix can be generated using the `binarize()` function.
#' @param parallel A logical value indicating whether to run the calculations in parallel. Defaults to FALSE.
#' @param lag.max The maximum lag to use in the cross-correlation function (CCF). Defaults to 1.
#' @param correlation_threshold The threshold value for filtering edges in the network (Pearsons coefficients go from -1 to +1). Set to "none" to disable filtering. Defaults to 0.5.
#' @return An `igraph` object representing the network of correlated cells.
#' @examples
#' binarized_calcium_matrix <- matrix(sample(c(0, 1), 1000, replace = TRUE), nrow = 10)
#' network <- make_network(binarized_calcium_matrix)
#' network_parallel <- make_network(binarized_calcium_matrix, parallel = TRUE, lag.max = 2, correlation_threshold = 0.4)
#' network_no_filter <- make_network(binarized_calcium_matrix, correlation_threshold = "none")
#' @export
#' @importFrom igraph graph.adjacency delete.edges E
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom utils txtProgressBar setTxtProgressBar
make_network <- function(binarized_calcium_matrix, parallel = FALSE, lag.max = 1, correlation_threshold = 0.5) {

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

  if (parallel) {
    # Set up parallel processing
    n_cores <- parallel::detectCores()
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("T.allcellst", "max_CCF", "lag.max"))
    parallel::clusterEvalQ(cl, {
      library(matrixStats)
    })

    # Function to compute max_CCF in parallel
    max_CCF_parallel <- function(i, j) {
      max_CCF(T.allcellst[, i], T.allcellst[, j], lag.max)
    }

    # Loop through upper triangular part of the matrix in parallel
    pb <- utils::txtProgressBar(min = 0, max = n_cols, style = 3)
    for (i in 1:(n_cols - 1)) {
      res <- parallel::parLapply(cl, (i + 1):n_cols, max_CCF_parallel, i = i)
      cmat.allcellst[i, (i + 1):n_cols] <- unlist(res)
      utils::setTxtProgressBar(pb, i)
    }

    # Clean up parallel processing
    parallel::stopCluster(cl)
    close(pb)
  } else {
    # Sequential computation
    pb <- utils::txtProgressBar(min = 0, max = n_cols, style = 3)
    for (i in 1:(n_cols - 1)) {
      for (j in (i + 1):n_cols) {
        cmat.allcellst[i, j] <- max_CCF(T.allcellst[, i], T.allcellst[, j], lag.max)
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  }

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
