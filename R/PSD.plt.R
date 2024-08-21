#' Plot Power Spectral Density (PSD)
#'
#' This function performs a Power Spectral Density (PSD) analysis on a calcium imaging matrix
#' and returns only the PSD plot. The user can specify whether to binarize the matrix before
#' the analysis and set the frame rate.
#'
#' @param calcium_matrix A matrix where each row represents a cell and each column represents a timepoint.
#' @param binarize A logical value indicating whether to binarize the calcium matrix before performing the PSD analysis. Defaults to TRUE.
#' @param frame_rate The frame rate of the calcium imaging data in frames per second. Defaults to 0.5 Hz (2 seconds per frame).
#' @return A ggplot object representing the PSD plot.
#' @examples
#' calcium_matrix <- matrix(runif(1000), nrow = 10)
#' PSD_plot <- PSD.plt(calcium_matrix, binarize = TRUE, frame_rate = 0.5)
#' print(PSD_plot)
#' @export
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs
#' @importFrom reshape2 melt
#' @importFrom dplyr summarise
PSD.plt <- function(calcium_matrix, binarize = TRUE, frame_rate = 0.5) {
  # Binarize the calcium matrix if requested
  if (binarize) {
    calcium_matrix <- binarize(calcium_matrix)
  }

  # Perform FFT on each cell's time series
  fft_result <- apply(calcium_matrix, 1, fft)

  # Calculate power spectral density (PSD) for each time series
  psd <- abs(fft_result)^2 / ncol(calcium_matrix) # Normalize by length of acquisition

  # Define frequency parameters
  fs <- frame_rate
  nyquist <- fs / 2
  lowest_freq <- 1 / ncol(calcium_matrix)
  freq <- seq(lowest_freq, nyquist, length.out = nrow(psd)) # Calculate frequency range

  # Convert PSD to data frame and melt for easier plotting
  psd.melt <- as.data.frame(psd)
  psd.melt$frequency <- freq
  psd.melt <- reshape2::melt(psd.melt, id.vars = "frequency", variable.name = "cell", value.name = "PSD")

  # Calculate the mean PSD for each frequency
  psd.mean <- dplyr::summarise(psd.melt, PSD_mean = mean(PSD), .by = "frequency")

  # Generate and return the PSD plot
  psd_plot <- ggplot2::ggplot(psd.mean, ggplot2::aes(x = frequency, y = PSD_mean)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Power Spectral Density (PSD)", x = "Frequency (Hz)", y = "PSD Mean")

  return(psd_plot)
}
