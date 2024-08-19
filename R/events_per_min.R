#' Calculate Event Frequency per Minute
#'
#' This function calculates the frequency of events per minute for each cell in a binarized calcium matrix.
#'
#' @param calcium_matrix_binary A binarized matrix where each row represents a cell and each column represents a timepoint. Can be created with binarize()
#' @param frame_rate The frame rate of the calcium imaging data (frames per second).
#' @param whole_population A logical value indicating whether to return the mean frequency of events per minute for all cells. Defaults to FALSE.
#' @return A numeric array representing the event frequency per minute for each cell, or a single numeric value representing the mean frequency if whole_population is TRUE.
#' @examples
#' binarized_data <- matrix(sample(c(0, 1), 100, replace = TRUE), nrow = 10)
#' frame_rate <- 30
#' event_frequency <- events_per_min(binarized_data, frame_rate)
#' mean_event_frequency <- events_per_min(binarized_data, frame_rate, TRUE)
#' @export
events_per_min <- function(calcium_matrix_binary, frame_rate, mean_all = FALSE) {
  # Count the number of events (spikes) for each cell
  event_counts <- rowSums(calcium_matrix_binary > 0)

  # Calculate the total recording time in minutes
  total_time_minutes <- ncol(calcium_matrix_binary) / frame_rate / 60

  # Calculate the event frequency per minute for each cell
  event_frequency <- event_counts / total_time_minutes

  # Return the mean frequency if whole_population is TRUE
  if (mean_all) {
    return(mean(event_frequency))
  } else {
    return(event_frequency)
  }
}
