
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Clustering

<!-- badges: start -->

[![R-CMD-check](https://github.com/ninjanick02/Clustering/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ninjanick02/Clustering/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ninjanick02/Clustering/graph/badge.svg)](https://app.codecov.io/gh/ninjanick02/Clustering)
<!-- badges: end -->

The **Clustering** package provides efficient implementations of the
k-means clustering algorithm using Lloyd's method. It includes both a
pure R implementation with vectorized operations and a high-performance
C++ implementation via Rcpp.

## Features

- **Two implementations**: Choose between pure R (`Kmeans()`) or
  optimized C++ (`Kmeans_cpp()`)
- **Vectorized operations**: The R implementation uses efficient matrix
  operations
- **Robust handling**: Gracefully handles empty clusters by
  re-initialization
- **Comprehensive output**: Returns cluster assignments, centers,
  within-cluster sum of squares, and iteration counts
- **Flexible initialization**: Accepts either the number of clusters or
  initial center coordinates
- **Utility functions**: Includes elbow method and silhouette score
  calculation

## Installation

You can install the development version of Clustering from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ninjanick02/Clustering")
```

## Example

Here's a basic example using the classic iris dataset:

``` r
library(Clustering)

# Prepare the data
X_iris <- as.matrix(iris[, 1:4])

# Run k-means clustering with k=3 clusters
set.seed(123)
result <- Kmeans(X_iris, centers = 3, max_iter = 100, tol = 1e-4)

# View the cluster centers
print(result$centers)

# View cluster assignments (first 10 observations)
print(result$cluster[1:10])

# Check convergence
cat("Converged in", result$iter, "iterations\n")
cat("Total within-cluster sum of squares:", result$tot.withinss, "\n")
```

### Using the C++ implementation

For larger datasets, the C++ implementation provides better performance:

``` r
# Run k-means with C++ implementation
set.seed(123)
result_cpp <- Kmeans_cpp(X_iris, centers = 3, max_iter = 100, tol = 1e-4)

# Results are in the same format
print(result_cpp$centers)
```

### Finding the optimal number of clusters

Use the elbow method to find the optimal k:

``` r
# Test different values of k
elbow_data <- elbow_method(X_iris, k_range = 2:6, n_init = 10)
print(elbow_data)

# Plot the elbow curve
plot(elbow_data$k, elbow_data$tot.withinss, 
     type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster SS",
     main = "Elbow Method for Optimal k")
```

### Evaluating clustering quality

Calculate silhouette scores to evaluate cluster quality:

``` r
# Perform clustering
result <- Kmeans(X_iris, centers = 3)

# Calculate silhouette scores
sil <- silhouette_score(X_iris, result$cluster)
cat("Average silhouette score:", sil$avg_silhouette, "\n")
```

## Performance

The package includes both R and C++ implementations to balance ease of
use and performance:

- **R implementation** (`Kmeans`): Uses vectorized matrix operations
  for good performance in pure R
- **C++ implementation** (`Kmeans_cpp`): Provides additional speed
  improvements for larger datasets

## Getting Help

If you encounter any issues or have questions:

- Report bugs at <https://github.com/ninjanick02/Clustering/issues>
- See the vignette for more detailed examples:
  `vignette("k_means_clustering", package = "Clustering")`
- Read the contributing guidelines in CONTRIBUTING.md
