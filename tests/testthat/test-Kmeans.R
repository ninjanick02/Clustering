library(testthat)
library(stats)
library(mclust)
library(Clustering)

test_that("R implementation matches stats::kmeans", {
  set.seed(123)
  X <- as.matrix(iris[, 1:4])
  k <- 3
  initial_indices <- sample.int(nrow(X), k)
  initial_centers_matrix <- X[initial_indices, , drop = FALSE]
  res_my_r <- Kmeans(X,
                     centers = initial_centers_matrix,
                     max_iter = 100,
                     tol = 1e-4)

  res_stats <- stats::kmeans(X,
                             centers = initial_centers_matrix,
                             iter.max = 100,
                             algorithm = "Lloyd")

  ari <- adjustedRandIndex(res_my_r$cluster, res_stats$cluster)
  expect_equal(ari, 1.0)#, tolerance = 1e-10)

  expect_equal(res_my_r$centers, res_stats$centers, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("C++ implementation matches stats::kmeans", {
  set.seed(123)
  X <- as.matrix(iris[, 1:4])
  k <- 3
  initial_indices <- sample.int(nrow(X), k)
  initial_centers_matrix <- X[initial_indices, , drop = FALSE]
  res_my_cpp <- Kmeans_cpp(X,
                           centers = initial_centers_matrix,
                           max_iter = 100,
                           tol = 1e-4)

  res_stats <- stats::kmeans(X,
                             centers = initial_centers_matrix,
                             iter.max = 100,
                             algorithm = "Lloyd")

  ari <- adjustedRandIndex(res_my_cpp$cluster, res_stats$cluster)
  expect_gt(ari, 0.99)#, tolerance = 1e-10)

  expect_equal(res_my_cpp$centers, res_stats$centers, tolerance = 1e-4, ignore_attr = TRUE)
})
