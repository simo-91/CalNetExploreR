#' Run Analysis Pipeline on Calcium Imaging Data
#'
#' This function runs a comprehensive analysis pipeline on calcium imaging data, including normalization, binarization,
#' population activity plotting, network creation and plotting, PCA analysis, power spectral density (PSD) analysis,
#' and degree distribution analysis. The function requires the user to provide coordinates for the cells.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param coordinates A data frame containing X and Y coordinates for each cell. Must include columns "X", "Y", and "Cell".
#' @param dendrogram A logical value indicating whether to include a dendrogram in the population activity plot. Defaults to FALSE.
#' @param correlation_threshold A numeric value specifying the threshold for filtering edges by weight in the network analysis. Set to "none" to disable filtering. Defaults to 0.3.
#' @param frame_rate A numeric value specifying the frame rate (in Hz) for the PSD analysis. Defaults to 0.5.
#' @return A list containing the results of each analysis step, including plots.
#' \itemize{
#'   \item \code{normalized_matrix}: The normalized calcium matrix.
#'   \item \code{binarized_matrix}: The binarized calcium matrix.
#'   \item \code{population_activity_plot}: The population activity plot.
#'   \item \code{network}: The network object created from the binarized matrix.
#'   \item \code{network_plot}: The plot of the network graph.
#'   \item \code{degree_plot}: The degree distribution plot of the network.
#'   \item \code{pca_plot}: The scree plot from PCA analysis.
#'   \item \code{psd_plot}: The power spectral density (PSD) plot.
#' }
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' coordinates <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10)
#' results <- pipeline(calcium_matrix, coordinates = coordinates)
#' @export

pipeline <- function(calcium_matrix, coordinates, dendrogram = FALSE, correlation_threshold = 0.3, frame_rate = 0.5) {

  # Ensure coordinates are provided and valid
  if (is.null(coordinates) || !all(c("X", "Y", "Cell") %in% colnames(coordinates))) {
    stop("Coordinates must be provided and must include columns 'X', 'Y', and 'Cell'.")
  }

  # Step 1: Normalize the calcium matrix
  normalized_matrix <- normalize(calcium_matrix)

  # Step 2: Binarize the normalized matrix
  binarized_matrix <- binarize(normalized_matrix)

  # Step 3: Population Activity Plotting
  pop_activity_plot <- population_activity.plt(binarized_matrix, binarize = FALSE, dendrogram = dendrogram)

  # Step 4: Network Creation
  network <- make_network(binarized_matrix, lag.max = 1, correlation_threshold = correlation_threshold)

  # Step 5: Network Plotting
  network_plot <- plot_network(graph = network, coordinates = coordinates, label = "communities", correlation_threshold = correlation_threshold)

  # Step 6: Degree Distribution Analysis
  degree_plot <- degrees(graph = network, plot = TRUE)

  # Step 7: PCA Analysis
  pca_result <- pca(normalized_matrix, binarize = FALSE, plot = TRUE)
  pca_plot <- pca_result

  # Step 8: PSD Analysis
  psd_plot <- PSD.plt(normalized_matrix, binarize = FALSE, frame_rate = frame_rate)

  # Print all plots individually
  print(pop_activity_plot)
  print(network_plot)
  print(degree_plot)
  print(pca_plot)
  print(psd_plot)

  # Return a list containing all results
  return(list(
    normalized_matrix = normalized_matrix,
    binarized_matrix = binarized_matrix,
    population_activity_plot = pop_activity_plot,
    network = network,
    network_plot = network_plot,
    degree_plot = degree_plot,
    pca_plot = pca_plot,
    psd_plot = psd_plot
  ))
}

