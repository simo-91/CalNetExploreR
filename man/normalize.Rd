% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize.R
\name{normalize}
\alias{normalize}
\title{Normalize Calcium Imaging Data}
\usage{
normalize(calcium_matrix)
}
\arguments{
\item{calcium_matrix}{A matrix where each row represents a cell and each column represents a timepoint.}
}
\value{
A normalized matrix where each cell's timeseries data is scaled to [0, 1].
}
\description{
This function normalizes the timeseries calcium data for each cell.
}
\examples{
data <- matrix(runif(100), nrow = 10)
normalized_data <- normalize(data)
}
