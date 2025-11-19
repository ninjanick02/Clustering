#' @useDynLib Clustering, .registration = TRUE
NULL

#' @importFrom Rcpp sourceCpp
NULL

#' Basic K-Means Clustering
#'
#' @description
#' Performs k-means clustering on a numeric matrix using Lloyd's algorithm.
#'
#' @param X A numeric matrix of data. Each row is an observation, and each
#'   column is a variable.
#' @param centers Either the number of clusters (k) or a matrix of initial
#'   cluster centers.
#' @param max_iter The maximum number of iterations allowed.
#' @param tol The tolerance for convergence. The algorithm stops if the sum of
#'   squared distances between old and new centers is less than this value.
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
#' @examples
#' # Use the iris dataset
#' X_iris <- as.matrix(iris[, 1:4])
#'
#' # Run k-means with k=3
#' result <- Kmeans(X_iris, centers = 3, max_iter = 100, tol = 1e-4)
#' print(result$centers)
#'
#' @export
Kmeans <- function(X, centers, max_iter = 100, tol = 1e-4) {

  # --- Input Validation ---
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

  # --- Initialization ---
  if (is.numeric(centers) && length(centers) == 1) {
    k <- as.integer(centers)
    if (k <= 0 || k > nrow(X)) {
      stop("k must be a positive integer less than or equal to the number of rows in X.")
    }
    initial_indices <- sample.int(nrow(X), k)
    centers_matrix <- X[initial_indices, , drop = FALSE]
  } else if (is.matrix(centers)) {
    centers_matrix <- centers
    k <- nrow(centers_matrix)
    if (ncol(X) != ncol(centers_matrix)) {
      stop("Dimensions of X and initial centers matrix do not match.")
    }
    if (anyNA(centers_matrix)) {
      stop("Initial centers matrix contains missing values.")
    }
  } else {
    stop("'centers' must be either an integer (k) or a matrix of initial centers.")
  }

  colnames(centers_matrix) <- colnames(X)
  cluster_assignments <- rep(0, nrow(X))

  # Iteration loop
  for (i in 1:max_iter) {

    # --- 1. Assignment Step (VECTORIZED) ---
    # This replaces your for-loop, which was the main bottleneck.
    # We use matrix algebra to calculate all squared distances at once.
    # The formula for squared Euclidean distance (a-b)^2 is (a^2 - 2ab + b^2)

    X_sq_sums <- rowSums(X^2)                     # a^2 (for each row in X)
    centers_sq_sums <- rowSums(centers_matrix^2) # b^2 (for each row in centers_matrix)
    prod_matrix <- X %*% t(centers_matrix)       # ab (the n x k matrix of dot products)

    dist_matrix <- -2 * prod_matrix
    dist_matrix <- sweep(dist_matrix, 2, centers_sq_sums, "+") # -2ab + b^2
    dist_matrix <- sweep(dist_matrix, 1, X_sq_sums, "+")     # -2ab + b^2 + a^2

    new_cluster_assignments <- apply(dist_matrix, 1, which.min)

    # --- 2. Convergence Check (Part 1) ---
    if (identical(new_cluster_assignments, cluster_assignments)) {
      break
    }
    cluster_assignments <- new_cluster_assignments

    # --- 3. Update Step ---
    # NOTE: This loop is still here. While it could be vectorized
    # (e.g., using rowsum), your loop-based approach is very
    # robust for handling empty clusters, which is good.
    # The real bottleneck was the Assignment Step, not this one.
    new_centers <- matrix(0, nrow = k, ncol = ncol(X))
    for (j in 1:k) {
      cluster_points <- X[cluster_assignments == j, , drop = FALSE]
      if (nrow(cluster_points) == 0) {
        warning(paste("Empty cluster in iteration", i, "- re-initializing cluster", j))
        new_centers[j, ] <- X[sample.int(nrow(X), 1), ]
      } else {
        new_centers[j, ] <- colMeans(cluster_points)
      }
    }

    # --- 4. Convergence Check (Part 2) ---
    center_change <- sum((new_centers - centers_matrix)^2)
    centers_matrix <- new_centers
    if (center_change < tol) {
      break
    }
  } # End of iteration loop

  # --- 5. Calculate WSS (VECTORIZED) ---
  # This replaces your final for-loop for better performance.

  # Get the center coordinates for each point's assigned cluster
  point_centers <- centers_matrix[cluster_assignments, ]

  # Calculate squared differences for all points at once
  point_sq_dists <- rowSums((X - point_centers)^2)

  # Sum the squared distances by cluster
  withinss <- tapply(point_sq_dists, cluster_assignments, sum)

  # Ensure all clusters are present, even if empty (with WSS = 0)
  withinss_full <- numeric(k)
  names(withinss_full) <- 1:k
  withinss_full[names(withinss)] <- withinss
  withinss_full[is.na(withinss_full)] <- 0

  tot_withinss <- sum(withinss_full)

  # --- 6. Return Results ---
  list(
    centers = centers_matrix,
    cluster = cluster_assignments,
    iter = i,
    withinss = withinss_full,
    tot.withinss = tot_withinss
  )
}
