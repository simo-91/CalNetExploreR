#' Extract the Cumulative Variance Explained by the Top 5 Principal Components
#'
#' This function extracts the cumulative variance explained by the top 5 principal components (PCs) from
#' a PCA result object. It returns the cumulative variance percentage up to the 5th principal component.
#'
#' @param pca_result A list object containing PCA results. It should include the eigenvalues of the PCA, specifically
#' the `eigenvalues` field with a `cumulative.variance.percent` component.
#'
#' @return A numeric value representing the cumulative variance explained by the first 5 principal components (PCs). If
#' the `pca_result` is \code{NULL} or does not contain eigenvalues, the function returns \code{NaN}.
#'
#' @examples
#' # Assuming `pca_result` is the output from a PCA analysis
#' # Example with a simulated pca_result
#' pca_result <- list(eigenvalues = list(cumulative.variance.percent = c(25, 40, 55, 65, 75, 80)))
#' get_top5pc_variance(pca_result)
#'
#' @export
get_top5pc_variance <- function(pca_result) {
  if (!is.null(pca_result$eigenvalues)) {
    # Extract cumulative variance for the top 5 PCs
    cumulative_variance <- pca_result$eigenvalues$cumulative.variance.percent
    # Return the cumulative variance for the 5th principal component
    return(cumulative_variance[5])
  } else {
    return(NaN)
  }
}
