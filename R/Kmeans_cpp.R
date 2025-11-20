#' @useDynLib Clustering, .registration = TRUE
NULL

#' K-Means Clustering (C++ Implementation)
#'
#' @description
#' Performs k-means clustering on a numeric matrix using Lloyd's algorithm
#' with a high-performance C++ backend via Rcpp. This implementation is
#' recommended for larger datasets where performance is critical.
#'
#' @param X A numeric matrix of data. Each row is an observation, and each
#'   column is a variable.
#' @param centers Either the number of clusters (k) or a matrix of initial
#'   cluster centers.
#' @param max_iter The maximum number of iterations allowed. Default is 100.
#' @param tol The tolerance for convergence. The algorithm stops if the sum of
#'   squared distances between old and new centers is less than this value.
#'   Default is 1e-4.
#'
#' @return
#' A list with the following components:
#' \item{centers}{A matrix of the final cluster centers.}
#' \item{cluster}{A vector of integers (from 1:k) indicating the cluster
#'   to which each point is allocated.}
#' \item{iter}{The number of iterations performed.}
#' \item{withinss}{A vector of the within-cluster sum of squares for each
#'   cluster.}
#' \item{tot.withinss}{The sum of all `withinss`.}
#'
#' @details
#' This function uses a C++ implementation of Lloyd's algorithm for improved
#' performance compared to the pure R implementation. The core iteration loop
#' runs in compiled C++ code, making it significantly faster for large datasets.
#'
#' Empty clusters are handled by re-initializing them to random data points.
#'
#' The algorithm stops when either:
#' - Cluster assignments no longer change, or
#' - The sum of squared distances between old and new centers is less than `tol`, or
#' - The maximum number of iterations (`max_iter`) is reached.
#'
#' @examples
#' # Use the iris dataset
#' X_iris = as.matrix(iris[, 1:4])
#'
#' # Run k-means with k=3
#' set.seed(123)
#' result = Kmeans_cpp(X_iris, centers = 3, max_iter = 100, tol = 1e-4)
#' print(result$centers)
#'
#' # Use custom initial centers
#' initial = X_iris[c(1, 50, 100), ]
#' result2 = Kmeans_cpp(X_iris, centers = initial)
#'
#' @seealso \code{\link{Kmeans}} for the pure R implementation
#' @export
Kmeans_cpp = function(X, centers, max_iter = 100, tol = 1e-4) {

  # --- Input Validation (R-side) ---
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }

  # Check for missing values
  if (anyNA(X)) {
    stop("X contains missing values (NA or NaN). Please remove or impute them before clustering.")
  }

  # Check for infinite values
  if (any(is.infinite(X))) {
    stop("X contains infinite values. Please handle them before clustering.")
  }

  if (nrow(X) < 1) {
    stop("X must have at least one row.")
  }

  if (ncol(X) < 1) {
    stop("X must have at least one column.")
  }

  if (max_iter <= 0) {
    stop("max_iter must be a positive integer.")
  }

  if (tol < 0) {
    stop("tol must be non-negative.")
  }

  # --- Initialization (R-side) ---
  if (is.numeric(centers) && length(centers) == 1) {
    # 'centers' is an integer k
    k = as.integer(centers)
    if (k <= 0 || k > nrow(X)) {
      stop("Number of clusters (k) must be a positive integer",
           " less than or equal to the number of rows in X.")
    }
    # Randomly select k rows from X to be the initial centers
    initial_indices = sample.int(nrow(X), k)
    initial_centers_matrix = X[initial_indices, , drop = FALSE]

  } else if (is.matrix(centers)) {
    # 'centers' is an initial centers matrix
    initial_centers_matrix = centers
    k = nrow(initial_centers_matrix)
    if (ncol(X) != ncol(initial_centers_matrix)) {
      stop("Dimensions of X and initial centers matrix do not match.")
    }
    if (anyNA(initial_centers_matrix)) {
      stop("Initial centers matrix contains missing values.")
    }

  } else {
    stop("'centers' must be either an integer (k) or a matrix of initial centers.")
  }

  # --- Call the C++ function for the heavy lifting ---
  # .Call() is not needed, Rcpp_Exports.R creates a nice R wrapper for us
  cpp_result = kmeans_cpp_loop(X, initial_centers_matrix, max_iter, tol)

  # --- Post-processing (R-side) ---
  # The C++ function returns the core components.
  # We calculate WSS here in R for simplicity.
  centers = cpp_result$centers
  cluster_assignments = cpp_result$cluster

  withinss = numeric(k)
  for (j in 1:k) {
    cluster_points = X[cluster_assignments == j, , drop = FALSE]
    if(nrow(cluster_points) > 0) {
      centered_points = sweep(cluster_points, 2, centers[j, ], "-")
      withinss[j] = sum(centered_points^2)
    }
  }

  tot_withinss = sum(withinss)

  # --- Return Results (same format as my_kmeans) ---
  list(
    centers = centers,
    cluster = cluster_assignments,
    iter = cpp_result$iter,
    withinss = withinss,
    tot.withinss = tot_withinss
  )
}
