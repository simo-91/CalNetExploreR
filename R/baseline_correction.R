#' Baseline correction using dF/F0
#'
#' Computes the baseline-corrected calcium signal using the dF/F0 method.
#' F0 is defined as the mean of the lowest quantile (default: 10 percent) of fluorescence values for each cell.
#'
#' @param calcium_matrix A numeric matrix (cells x timepoints).
#' @param quantile_threshold Numeric value between 0 and 1 indicating which quantile to use as F0 (default: 0.1 for lowest 10 percent).
#'
#' @return A matrix of the same dimensions, containing dF/F0 values.
#'
#' @examples
#' mat <- matrix(runif(50, 100, 200), nrow = 5)
#' baseline_correction(mat)
#'
#' @export
baseline_correction <- function(calcium_matrix, quantile_threshold = 0.1) {
  corrected <- t(apply(calcium_matrix, 1, function(x) {
    f0 <- mean(x[x <= stats::quantile(x, quantile_threshold, na.rm = TRUE)], na.rm = TRUE)
    if (is.na(f0) || f0 == 0) {
      return(rep(0, length(x)))  # fallback to zeros if baseline is invalid
    }
    (x - f0) / f0
  }))
  rownames(corrected) <- rownames(calcium_matrix)
  colnames(corrected) <- colnames(calcium_matrix)
  corrected
}
