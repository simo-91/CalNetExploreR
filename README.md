# CalNetExploreR

**CalNetExploreR** is an R package designed for comprehensive calcium imaging data analysis, including normalization, binarization, network creation, and various advanced visualizations. The package provides an integrated pipeline to streamline the analysis of calcium imaging data from live imaging experiments.

## Installation

You can install the development version of CalNetExploreR from GitHub using the following commands:

Install devtools if you haven't already:

`install.packages("devtools")`

The `CalNetExploreR` package relies on several R packages for its functionality. Make sure to install these dependencies before using `CalNetExploreR`. You can install these packages using the command `install.packages("package_name")`.

Here is the list of dependencies:

- **ggplot2**: For data visualization and plotting graphs.
- **reshape2**: For data reshaping and melting data frames.
- **ggdendro**: For creating dendrogram plots.
- **cowplot**: For combining multiple plots into a single plot.
- **ggpubr**: For publication-ready plots with additional theme options.
- **grid**: For viewport manipulation and custom plot arrangements.
- **igraph**: For creating and analyzing networks and graphs.
- **ggraph**: For plotting networks and graphs with ggplot2.
- **factoextra**: For visualizing PCA results and extracting eigenvalues.
- **dplyr**: For data manipulation and summarization.
- **RColorBrewer**: For color palettes used in plots.

### Installation of Dependencies

To ensure all dependencies are installed, you can run the following command:

`install.packages(c("ggplot2", "reshape2", "ggdendro", "cowplot", "ggpubr", 
                   "grid", "ggraph", "igraph", "factoextra", "dplyr", "RColorBrewer"))`

## Install CalNetExploreR from GitHub

`devtools::install_github("simo-91/CalNetExploreR")`

## Usage

To use the main functions you need a matrix containing your timeseries data. Every row should be a different cell/timeseries and it should looks something like this:

| Timepoint n~0~ | Timepoint 2 | Timepoint 3 | ... | Timepoint n~max~ |
|:--------------:|:-----------:|:-----------:|:---:|:----------------:|
|      0.12      |    0.23     |    0.35     | ... |       0.31       |
|      0.11      |    0.19     |    0.22     | ... |       0.40       |
|      0.09      |    0.15     |    0.18     | ... |       0.42       |
|      0.14      |    0.21     |    0.30     | ... |       0.37       |
|      0.13      |    0.20     |    0.25     | ... |       0.43       |
|      0.08      |    0.12     |    0.16     | ... |       0.39       |
|      0.10      |    0.18     |    0.23     | ... |       0.40       |

Some of the functions require cell coordinates (X and Y) stored in a matrix. This should have X, Y and Cell ID columns, e.g.:

|  X   |  Y   | Cell ID |
|:----:|:----:|:-------:|
| 12.5 | 23.8 |    1    |
| 14.1 | 25.3 |    2    |
| 10.2 | 22.7 |    3    |
| 13.7 | 20.9 |    4    |
| 15.4 | 18.5 |    5    |
| 11.8 | 19.3 |    6    |
| 16.2 | 21.4 |    7    |
| ...  | ...  |   ...   |

## Usage example:

```         
# Load the package

library(CalNetExploreR)

# Example data: Generate a random calcium matrix and coordinates
    
calcium_matrix <- matrix(runif(1000), nrow = 10)  
coordinates <- data.frame(X = runif(10), Y = runif(10), Cell = 1:10)

# Population Activity Plot with a dendrogram
population_plot <- population_activity.plt(binarized_calcium_matrix = results$binarized_matrix, binarize = FALSE, dendrogram = TRUE)
print(population_plot)

# PCA Analysis with a scree plot
pca_result <- pca(calcium_matrix = results$normalized_matrix, binarize = FALSE, plot = TRUE)
print(pca_result)

# Create a network from the binarized matrix
network <- make_network(binarized_calcium_matrix = results$binarized_matrix, lag.max = 1, correlation_threshold = 0.3)

# Plot the network with community labels
network_plot <- plot_network(graph = network, coordinates = coordinates, label = "communities", correlation_threshold = 0.3)
print(network_plot)

# Degree Distribution Plot
degree_plot <- degrees(graph = network, plot = TRUE)
print(degree_plot)

# Power Spectral Density Analysis Plot
psd_plot <- PSD.plt(calcium_matrix = results$normalized_matrix, binarize = FALSE, frame_rate = 0.5)
print(psd_plot)

# Calculate events per minute
events_per_min_results <- events_per_min(calcium_matrix = results$binarized_matrix)
print(events_per_min_results)
```

## Main Features

CalNetExploreR offers the following functionalities:

-   Normalization: Standardize calcium imaging data to a common scale.
-   Binarization: Convert normalized data to binary states for activity analysis.
-   Population Activity Plotting: Hierarchical clustering to sort and display cell activity over time with optional dendrogram clustering + percentage of coactive cells
-   Network Analysis: Create and visualize a network of cell interactions based on cross-correlation and lag.
-   Principal Component Analysis (PCA): Analyze data variance and display scree plots of eigenvalues.
-   Power Spectral Density (PSD) Analysis: Perform frequency analysis on calcium signals. Degree Distribution Analysis: Analyze the degree distribution of nodes within the network.

For more information, please refer to the manual (CalNetExploreR.pdf)

## License

This package is open-source and available under the MIT License.
