


#' Main Function: k-means clustering for count data
#'
#' @export
kmeansPoi <- function(data, k, min_mu = 1e-4, verbose = FALSE, debug = FALSE){

  # Parameter validation
  Y <- data


  # Initialize
  size_factors <- matrixStats::colMeans2(Y)
  # centers <- Y[, sample(seq_len(ncol(Y)), size = k)] + min_mu
  # centers <- matrix(NA, ncol = k, nrow = nrow(Y))
  # random_assign <- sample(seq_len(ncol(Y)), size = k * 100, replace = FALSE)
  # for(ki in seq_len(k)){
  #   sel <- random_assign[seq_along(random_assign) %% k == ki -1]
  #   centers[,ki] <- pmax(rowSums2(Y, cols = sel) / sum(size_factors[sel]), min_mu)
  # }
  centers <- kmeans_pp_initialization(Y, size_factors, k = k, min_mu = min_mu)

  # Call clustering
  if(debug){
    run_kmeans_r(Y, size_factors, centers,
               min_mu = min_mu, max_iter = 100,
               tolerance = 1e-8, verbose = verbose)
  }else{
    run_kmeans(Y, size_factors, centers,
               min_mu = min_mu, max_iter = 100,
               tolerance = 1e-8, verbose = verbose)
  }

}


kmeansPoi2 <- function(data, k, min_mu = 1e-4, verbose = FALSE, debug = FALSE){

  # Parameter validation
  Y <- data


  # Initialize
  size_factors <- matrixStats::colMeans2(Y)
  # centers <- Y[, sample(seq_len(ncol(Y)), size = k)] + min_mu
  # centers <- matrix(NA, ncol = k, nrow = nrow(Y))
  # random_assign <- sample(seq_len(ncol(Y)), size = k * 100, replace = FALSE)
  # for(ki in seq_len(k)){
  #   sel <- random_assign[seq_along(random_assign) %% k == ki -1]
  #   centers[,ki] <- pmax(rowSums2(Y, cols = sel) / sum(size_factors[sel]), min_mu)
  # }
  centers <- kmeans_pp_initialization2(Y, size_factors, k = k, min_mu = min_mu)

  # Call clustering
  if(debug){
    run_kmeans_r(Y, size_factors, centers,
                 min_mu = min_mu, max_iter = 100,
                 tolerance = 1e-8, verbose = verbose)
  }else{
    run_kmeans2(Y, size_factors, centers,
               min_mu = min_mu, max_iter = 100,
               tolerance = 1e-8, verbose = verbose)
  }

}
