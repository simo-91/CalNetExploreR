#' Perform PCA on Calcium Imaging Data
#'
#' This function performs Principal Component Analysis (PCA) on the calcium imaging data and optionally plots the scree plot of the eigenvalues.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param binarize A logical value indicating whether to binarize the calcium matrix before performing PCA. Defaults to TRUE.
#' @param plot A logical value indicating whether to plot the scree plot of eigenvalues. Defaults to TRUE.
#' @return If plot = FALSE, returns a list containing the PCA results and eigenvalues. If plot = TRUE, displays the scree plot and does not return anything.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' pca_results <- pca(calcium_matrix, binarize = TRUE, plot = TRUE)
#' @export
#' @importFrom stats prcomp
#' @importFrom factoextra fviz_eig get_eigenvalue
#' @importFrom ggpubr ggpar
#' @import ggplot2
pca <- function(calcium_matrix, binarize = TRUE, plot = TRUE) {
  # Binarize the calcium matrix if requested
  if (binarize) {
    calcium_matrix <- binarize(normalize(calcium_matrix))
  }

  # Perform PCA using the stats package
  pca_result <- stats::prcomp(calcium_matrix, center = TRUE, scale. = FALSE)

  # Extract eigenvalues
  pca_eigenvalues <- factoextra::get_eigenvalue(pca_result)

  # Plot the scree plot if requested
  if (plot) {
    scree_plot <- factoextra::fviz_eig(pca_result, ncp = 100)
    scree_plot <- ggpubr::ggpar(scree_plot, title = ggplot2::element_blank())
    print(scree_plot)  # Explicitly print the plot
  } else {
    # Return a list containing the PCA result and eigenvalues
    return(list(pca_result = pca_result, eigenvalues = pca_eigenvalues))
  }
}
