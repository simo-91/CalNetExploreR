#' Run Analysis Pipeline on Calcium Imaging Data
#'
#' This function runs a comprehensive analysis pipeline on calcium imaging data, including normalization, binarization,
#' population activity plotting, network creation and plotting, PCA analysis, power spectral density (PSD) analysis,
#' degree distribution analysis, and various network metrics calculations (e.g., clustering coefficient, global efficiency, and event frequency).
#'
#' In the network plotting step, if `cell_labels = TRUE`, nodes are labeled with the community numbers they belong to by setting
#' `cell_ID = "communities"` in the `plot_network()` function. This facilitates visualization of community structures within the network.
#'
#' If `cell_labels = FALSE`, the function skips all label-based calculations, including the check for the `Label` column in the `coordinates`
#' data frame and the computation of labeled-to-unlabeled connections.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param coordinates A data frame containing spatial coordinates for each cell. Must include columns "X", "Y", and "Cell".
#'        If `cell_labels = TRUE`, must also include a "Label" column.
#' @param dendrogram Logical. Whether to include a dendrogram in the population activity plot. Defaults to `FALSE`.
#' @param correlation_threshold Numeric. Threshold for filtering edges by weight in the network analysis.
#'        Set to `"none"` to disable filtering. Defaults to `0.3`.
#' @param frame_rate Numeric. The frame rate (in Hz) for the PSD analysis and event frequency calculation. Defaults to `0.5`.
#' @param lag.max Numeric. The maximum lag used in the network creation step. Defaults to `1`.
#' @param big_community_min_members Integer. Minimum number of members required for a community to be considered "big". Defaults to `5`.
#' @param samplename Character string. The name of the sample, used to name saved plot image files.
#' @param clustering_method Character string. The hierarchical clustering method to be used in the population activity plot.
#'        Options: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid". Defaults to `"ward.D2"`.
#' @param cell_labels Logical. Whether the `coordinates` data frame contains a "Label" column and whether label-based
#'        analysis (like labeled-to-unlabeled connection metrics) should be performed. Defaults to `TRUE`.
#'
#' @return A named list containing all intermediate outputs and calculated metrics, including:
#' \describe{
#'   \item{normalized_matrix}{Normalized calcium activity matrix}
#'   \item{binarized_matrix}{Binarized activity matrix}
#'   \item{population_activity_plot}{ggplot object of population raster plot}
#'   \item{network}{igraph object of the constructed network}
#'   \item{network_plot}{ggplot object of the network plot}
#'   \item{degree_plot}{ggplot object showing degree distribution}
#'   \item{pca_plot}{ggplot object of PCA}
#'   \item{psd_plot}{ggplot object of power spectral density plot}
#'   \item{top5pc_variance}{Variance explained by the top 5 principal components}
#'   \item{clustering_coefficient}{Network transitivity (clustering coefficient)}
#'   \item{global_efficiency}{Global efficiency of the network}
#'   \item{total_connections_labeled_to_unlabeled}{Number of edges from labeled to unlabeled cells (if `cell_labels = TRUE`)}
#'   \item{total_possible_connections_labeled_to_unlabeled}{Possible such edges (denominator)}
#'   \item{proportion_labeled_to_unlabeled}{Proportion of labeled-to-unlabeled connections}
#'   \item{total_connections_labeled_to_labeled}{Edges between labeled cells}
#'   \item{total_possible_connections_labeled_to_labeled}{Possible labeled-labeled connections}
#'   \item{proportion_labeled_to_labeled}{Proportion of labeled-labeled connections}
#'   \item{event_frequency_per_min}{Vector of per-cell event frequencies (events/min)}
#'   \item{mean_event_frequency_per_min}{Average event frequency (events/min)}
#' }
#'
#' @export
#' @importFrom ggplot2 ggsave

pipeline <- function(calcium_matrix, coordinates, dendrogram = FALSE, correlation_threshold = 0.3,
                     frame_rate = 0.5, lag.max = 1, big_community_min_members = 5,
                     samplename = "sample", clustering_method = "ward.D2", cell_labels = TRUE) {

  # Ensure coordinates are provided and valid
  required_columns <- c("X", "Y", "Cell")
  if (cell_labels) {
    required_columns <- c(required_columns, "Label")
  }

  if (is.null(coordinates) || !all(required_columns %in% colnames(coordinates))) {
    stop(paste("Coordinates must include columns:", paste(required_columns, collapse = ", ")))
  }

  # Step 1: Normalize the calcium matrix
  normalized_matrix <- normalize(calcium_matrix)

  # Step 2: Binarize the normalized matrix
  binarized_matrix <- binarize(normalized_matrix)

  # Step 3: Population Activity Plotting
  pop_activity_plot <- population_activity(binarized_matrix,
                                           binarize = TRUE,
                                           dendrogram = dendrogram)

  # Step 4: Network Creation
  network <- make_network(binarized_matrix, lag.max = lag.max, correlation_threshold = correlation_threshold)

  # Step 5: Network Plotting
  # If cell_labels = TRUE, label by communities; otherwise, no labels
  network_plot <- plot_network(
    graph = network,
    coordinates = coordinates,
    label = "communities",
    cell_ID = "none",
    correlation_threshold = correlation_threshold
  )

  # Step 6: Degree Distribution Analysis
  degree_plot <- degrees(graph = network, plot = TRUE)

  # Step 7: PCA Analysis
  pca_result <- pca(binarized_matrix, binarize = FALSE, plot = TRUE)
  pca_plot <- pca_result
  pca_result <- pca(binarized_matrix, binarize = FALSE, plot = FALSE)

  # Step 8: PSD Analysis
  psd_plot <- PSD.plt(binarized_matrix, binarize = FALSE, frame_rate = frame_rate)

  # Step 9: Calculate variance explained by top 5 PCs
  top5pc_variance <- get_top5pc_variance(pca_result)

  # Step 10: Clustering Coefficient and Global Efficiency
  clustering_coefficient <- igraph::transitivity(network)
  global_efficiency <- igraph::global_efficiency(network)

  # Step 11: Calculate Labeled-to-Unlabeled Connections
  if (cell_labels) {
    correlation_matrix <- igraph::as_adjacency_matrix(network, attr = "weight", sparse = FALSE)
    subset_conn_results <- subset_connections(correlation_matrix, coordinates, correlation_threshold = correlation_threshold)

    total_connections_labeled_to_unlabeled <- subset_conn_results$total_connections_labeled_to_unlabeled
    total_possible_connections_labeled_to_unlabeled <- subset_conn_results$total_possible_connections_labeled_to_unlabeled
    proportion_labeled_to_unlabeled <- subset_conn_results$proportion_labeled_to_unlabeled

    total_connections_labeled_to_labeled <- subset_conn_results$total_connections_labeled_to_labeled
    total_possible_connections_labeled_to_labeled <- subset_conn_results$total_possible_connections_labeled_to_labeled
    proportion_labeled_to_labeled <- subset_conn_results$proportion_labeled_to_labeled
  } else {
    total_connections_labeled_to_unlabeled <- NA
    total_possible_connections_labeled_to_unlabeled <- NA
    proportion_labeled_to_unlabeled <- NA
    total_connections_labeled_to_labeled <- NA
    total_possible_connections_labeled_to_labeled <- NA
    proportion_labeled_to_labeled <- NA
  }

  # Step 12: Event Frequency
  event_frequency_per_min <- events_per_min(binarized_matrix, frame_rate)
  mean_event_frequency_per_min <- events_per_min(binarized_matrix, frame_rate, mean_all = TRUE)

  # Save plots
  ggplot2::ggsave(paste0(samplename, "_population_activity_plot.png"), plot = pop_activity_plot)
  ggplot2::ggsave(paste0(samplename, "_network_plot.png"), plot = network_plot)
  ggplot2::ggsave(paste0(samplename, "_degree_plot.png"), plot = degree_plot)
  ggplot2::ggsave(paste0(samplename, "_pca_plot.png"), plot = pca_plot)
  ggplot2::ggsave(paste0(samplename, "_psd_plot.png"), plot = psd_plot)

  # Return all results
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
