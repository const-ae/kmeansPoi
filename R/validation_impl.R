poisson_deviance_r <- function(y, mu){
  ifelse(y == 0, 2 * mu, pmax(2 * (y * log(y / mu) - y + mu), 0))
}


kmeans_pp_initialization_r <- function(Y, size_factors, k, min_mu){
  cl_center_indices <- rep(NA, k)
  sel_c1 <- sample(seq_len(ncol(Y)), 1)
  cl_center_indices[1] <- sel_c1

  clusters <- matrix(NA, nrow = nrow(Y), ncol = k)
  clusters[,1] <- pmax(rowSums2(Y, cols = sel_c1) / size_factors[sel_c1], min_mu)
  dev <- apply(Y, 2, function(y) sum(poisson_deviance_r(y, mu =  clusters[,1])))

  for(ki in seq(2, k)){
    weights <- dev^2 / sum(dev^2)
    indices <- setdiff(seq_len(ncol(Y)), cl_center_indices)
    sel_ci <- sample(indices, prob = weights[indices], 1)
    cl_center_indices[ki] <- sel_ci
    clusters[,ki] <- pmax(rowSums2(Y, cols = sel_ci) / size_factors[sel_ci], min_mu)

    if(ki == k){
      break
    }
    new_dev <- apply(Y, 2, function(y) sum(poisson_deviance_r(y, mu =  clusters[,ki])))
    dev <- pmin(dev, new_dev)
  }

  clusters
}



run_kmeans_r <- function(Y, size_factors, centers,
                         min_mu, max_iter, tolerance, verbose){

  n_samples <- ncol(Y)
  k <- ncol(centers)

  Dev <- matrix(NA, nrow = n_samples, ncol = k)
  old_centers <- centers

  err_lr <- Inf
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
    # err <- sum((centers - old_centers)^2)
    # if(verbose){
    #   message("iter: ", iter, "\terror: ",  err)
    # }
    # if(err < tolerance){
    #   if(verbose){
    #     message("Converged, after ", iter, " iterations. Error: ", err)
    #   }
    #   break
    # }
    err <- sum(sapply(seq_len(ncol(Y)), function(idx) Dev[idx, assign[idx]]))
    if(verbose){
      message("iter: ", iter, "\terror: ",  err)
    }
    if(abs(err - err_lr) / (err + 0.1) < 1e-5){
      message("Converged, after ", iter, " iterations. Error: ", err)
      break
    }
    err_lr <- err
    old_centers <- centers
  }

  list(centers = centers, cluster = assign, iterations = iter)
}
