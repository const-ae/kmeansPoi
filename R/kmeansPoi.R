


#' Main Function: k-means clustering for count data
#'
#' @export
kmeansPoi <- function(data, k, min_mu = 1e-4, verbose = FALSE){

  # Parameter validation
  Y <- data


  # Initialize
  size_factors <- matrixStats::colMeans2(Y)
  centers <- Y[, sample(seq_len(ncol(Y)), size = k)] + min_mu


  # Call clustering
  run_kmeans(Y, size_factors, centers,
             min_mu = min_mu, max_iter = 100,
             tolerance = 1e-8, verbose = verbose)

}
