#include <Rcpp.h>
using namespace Rcpp;

// Helper function to calculate squared Euclidean distances
// and find the closest cluster center.
// (We generally don't export helpers, so no tag needed here,
// unless you want to test this function specifically in R).
IntegerVector find_closest_centers(const NumericMatrix& X, const NumericMatrix& centers) {
  int n_obs = X.nrow();
  int k = centers.nrow();
  int n_vars = X.ncol();
  IntegerVector assignments(n_obs);

  for (int i = 0; i < n_obs; ++i) {
    double min_dist = R_PosInf;
    int best_cluster = 0;

    for (int j = 0; j < k; ++j) {
      double dist = 0.0;
      for (int l = 0; l < n_vars; ++l) {
        double diff = X(i, l) - centers(j, l);
        dist += diff * diff;
      }

      if (dist < min_dist) {
        min_dist = dist;
        best_cluster = j + 1; // 1-based indexing for R
      }
    }
    assignments[i] = best_cluster;
  }

  return assignments;
}

// [[Rcpp::export]]
List kmeans_cpp_loop(const NumericMatrix& X, NumericMatrix centers, int max_iter, double tol) {
  int n_obs = X.nrow();
  int k = centers.nrow();
  int n_vars = X.ncol();

  NumericMatrix new_centers(k, n_vars);
  IntegerVector assignments(n_obs);
  IntegerVector old_assignments(n_obs);

  for (int i = 0; i < max_iter; ++i) {

    // 1. Assignment Step
    assignments = find_closest_centers(X, centers);

    // Check for convergence (if assignments didn't change)
    if (is_true(all(assignments == old_assignments))) {
      return List::create(
        _["centers"] = centers,
        _["cluster"] = assignments,
        _["iter"] = i
      );
    }
    old_assignments = clone(assignments);

    // 2. Update Step
    // Reset new_centers and cluster counts
    new_centers.fill(0.0);
    IntegerVector counts(k);

    for (int j = 0; j < n_obs; ++j) {
      int cluster_idx = assignments[j] - 1; // 0-based index
      counts[cluster_idx]++;
      for (int l = 0; l < n_vars; ++l) {
        new_centers(cluster_idx, l) += X(j, l);
      }
    }

    // Divide by counts to get the mean
    for (int j = 0; j < k; ++j) {
      if (counts[j] > 0) {
        for (int l = 0; l < n_vars; ++l) {
          new_centers(j, l) /= counts[j];
        }
      } else {
        // Handle empty cluster: re-initialize to a random point
        int rand_idx = floor(R::runif(0, n_obs));
        for(int l = 0; l < n_vars; ++l) {
          new_centers(j, l) = X(rand_idx, l);
        }
      }
    }

    // Check for tolerance-based convergence
    double center_change = 0.0;
    for (int j = 0; j < k; ++j) {
      for (int l = 0; l < n_vars; ++l) {
        double diff = new_centers(j, l) - centers(j, l);
        center_change += diff * diff;
      }
    }

    centers = clone(new_centers); // Update centers

    if (center_change < tol) {
      return List::create(
        _["centers"] = centers,
        _["cluster"] = assignments,
        _["iter"] = i + 1
      );
    }
  }

  // Reached max_iter
  return List::create(
    _["centers"] = centers,
    _["cluster"] = assignments,
    _["iter"] = max_iter
  );
}
