#' Run Analysis Pipeline on Calcium Imaging Data
#'
#' This function runs a comprehensive analysis pipeline on calcium imaging data, including normalization, binarization,
#' population activity plotting, network creation and plotting, PCA analysis, power spectral density (PSD) analysis,
#' degree distribution analysis, and various network metrics calculations (e.g., clustering coefficient, global efficiency, and event frequency).
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param coordinates A data frame containing X and Y coordinates for each cell. Must include columns "X", "Y", "Cell", and "Label".
#' @param dendrogram A logical value indicating whether to include a dendrogram in the population activity plot. Defaults to FALSE.
#' @param correlation_threshold A numeric value specifying the threshold for filtering edges by weight in the network analysis. Set to "none" to disable filtering. Defaults to 0.3.
#' @param frame_rate A numeric value specifying the frame rate (in Hz) for the PSD analysis and event frequency calculation. Defaults to 0.5.
#' @param lag.max A numeric value specifying the lag to be used in the network creation step. Defaults to 1.
#' @param big_community_min_members An integer specifying the minimum number of members for a community to be considered "big." Defaults to 5.
#' @param samplename A character string specifying the name of the sample. Used to name saved plot images.
#' @return A list containing the results of each analysis step, including plots and calculated features.
#' @export
#' @importFrom ggplot2 ggsave
pipeline <- function(calcium_matrix, coordinates, dendrogram = FALSE, correlation_threshold = 0.3, frame_rate = 0.5, lag.max = 1, big_community_min_members = 5, samplename = "sample") {

  # Ensure coordinates are provided and valid
  if (is.null(coordinates) || !all(c("X", "Y", "Cell", "Label") %in% colnames(coordinates))) {
    stop("Coordinates must be provided and must include columns 'X', 'Y', 'Cell', and 'Label'.")
  }

  # Step 1: Normalize the calcium matrix
  normalized_matrix <- normalize(calcium_matrix)

  # Step 2: Binarize the normalized matrix
  binarized_matrix <- binarize(normalized_matrix)

  # Step 3: Population Activity Plotting
  pop_activity_plot <- population_activity.plt(binarized_matrix, binarize = FALSE, dendrogram = dendrogram)

  # Step 4: Network Creation
  network <- make_network(binarized_matrix, lag.max = lag.max, correlation_threshold = correlation_threshold)

  # Step 5: Network Plotting
  network_plot <- plot_network(graph = network, coordinates = coordinates, label = "communities", correlation_threshold = correlation_threshold)

  # Step 6: Degree Distribution Analysis
  degree_plot <- degrees(graph = network, plot = TRUE)

  # Step 7: PCA Analysis (applies to binarized matrix)
  pca_result <- pca(binarized_matrix, binarize = FALSE, plot = TRUE)
  pca_plot <- pca_result
  pca_result <- pca(binarized_matrix, binarize = FALSE, plot = FALSE)

  # Step 8: PSD Analysis (applies to binarized matrix)
  psd_plot <- PSD.plt(binarized_matrix, binarize = FALSE, frame_rate = frame_rate)

  # Step 9: Calculate variance explained by the top 5 principal components
  top5pc_variance <- get_top5pc_variance(pca_result)

  # Step 10: Calculate Clustering Coefficient and Global Efficiency
  clustering_coefficient <- transitivity(network)
  global_efficiency <- global_efficiency(network)

  # Step 11: Calculate Labeled-to-NonLabeled Connections using subset_connections()
  # First, obtain the correlation matrix used in network creation
  correlation_matrix <- as_adjacency_matrix(network, attr = "weight", sparse = FALSE)

  # Use the subset_connections function
  subset_conn_results <- subset_connections(correlation_matrix, coordinates, correlation_threshold = correlation_threshold)

  # Extract the required metrics
  total_connections_labeled_to_unlabeled <- subset_conn_results$total_connections_labeled_to_unlabeled
  total_possible_connections_labeled_to_unlabeled <- subset_conn_results$total_possible_connections_labeled_to_unlabeled
  proportion_labeled_to_unlabeled <- subset_conn_results$proportion_labeled_to_unlabeled

  total_connections_labeled_to_labeled <- subset_conn_results$total_connections_labeled_to_labeled
  total_possible_connections_labeled_to_labeled <- subset_conn_results$total_possible_connections_labeled_to_labeled
  proportion_labeled_to_labeled <- subset_conn_results$proportion_labeled_to_labeled

  # Step 12: Calculate Event Frequency (events per minute)
  event_frequency_per_min <- events_per_min(binarized_matrix, frame_rate)
  mean_event_frequency_per_min <- events_per_min(binarized_matrix, frame_rate, mean_all = TRUE)

  # Save plots as images using samplename
  ggplot2::ggsave(paste0(samplename, "_population_activity_plot.png"), plot = pop_activity_plot)
  ggplot2::ggsave(paste0(samplename, "_network_plot.png"), plot = network_plot)
  ggplot2::ggsave(paste0(samplename, "_degree_plot.png"), plot = degree_plot)
  ggplot2::ggsave(paste0(samplename, "_pca_plot.png"), plot = pca_plot)
  ggplot2::ggsave(paste0(samplename, "_psd_plot.png"), plot = psd_plot)

  # Return a list containing all results and features
  return(list(
    normalized_matrix = normalized_matrix,
    binarized_matrix = binarized_matrix,
    population_activity_plot = pop_activity_plot,
    network = network,
    network_plot = network_plot,
    degree_plot = degree_plot,
    pca_plot = pca_plot,
    psd_plot = psd_plot,
    top5pc_variance = top5pc_variance,
    clustering_coefficient = clustering_coefficient,
    global_efficiency = global_efficiency,
    total_connections_labeled_to_unlabeled = total_connections_labeled_to_unlabeled,
    total_possible_connections_labeled_to_unlabeled = total_possible_connections_labeled_to_unlabeled,
    proportion_labeled_to_unlabeled = proportion_labeled_to_unlabeled,
    total_connections_labeled_to_labeled = total_connections_labeled_to_labeled,
    total_possible_connections_labeled_to_labeled = total_possible_connections_labeled_to_labeled,
    proportion_labeled_to_labeled = proportion_labeled_to_labeled,
    event_frequency_per_min = event_frequency_per_min,
    mean_event_frequency_per_min = mean_event_frequency_per_min
  ))
}
