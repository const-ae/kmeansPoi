poisson_deviance_r <- function(y, mu){
  ifelse(y == 0, 2 * mu, pmax(2 * (y * log(y / mu) - y + mu), 0))
}


run_kmeans_r <- function(Y, size_factors, centers,
                         min_mu, max_iter, tolerance, verbose){

  n_samples <- ncol(Y)
  k <- ncol(centers)

  Dev <- matrix(NA, nrow = n_samples, ncol = k)
  old_centers <- centers

  for(iter in seq_len(max_iter)){

    # Assignments
    for(ki in seq_len(k)){
      for(si in seq_len(n_samples)){
        Dev[si, ki] <- sum(poisson_deviance_r(Y[,si], mu = centers[,ki] * size_factors[si]))
      }
    }

    assign <- apply(Dev, 1, which.min)
    # print(table(assign))

    # Update Center
    for(ki in seq_len(k)){
      sel <- assign == ki
      centers[,ki] <- pmax(rowSums2(Y, cols = sel) / sum(size_factors[sel]), min_mu)
    }
    err <- sum((centers - old_centers)^2)
    if(verbose){
      message("iter: ", iter, "\terror: ",  err)
    }
    if(err < tolerance){
      if(verbose){
        message("Converged, after ", iter, " iterations. Error: ", err)
      }
      break
    }
    err_lr <- err
    old_centers <- centers
  }

  list(centers = centers, cluster = assign, iterations = iter)
}
