% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binarize.R
\name{binarize}
\alias{binarize}
\title{Binarize Calcium Imaging Data}
\usage{
binarize(
  calcium_matrix,
  cutoff_func = function(x) {
     th <- 2 * sd(x)
     x[x <= th] <- 0
     x[x > th] <-
    1
     return(x)
 }
)
}
\arguments{
\item{calcium_matrix}{A matrix where each row represents a cell and each column represents a timepoint.}

\item{cutoff_func}{A function to determine the threshold for binarizing the data. Default is twice the standard deviation of each cell.}
}
\value{
A binary matrix where each cell's timeseries data is converted to 0 or 1 based on the cutoff function.
}
\description{
This function binarizes the timeseries calcium data for each cell using a specified cutoff function.
}
\examples{
data <- matrix(runif(100), nrow = 10)
binary_data <- binarize(data)
}
