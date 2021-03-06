#include <Rcpp.h>
using namespace Rcpp;






inline double poisson_deviance_impl(const double y, const double mu){
  if(y == 0){
    return 2.0 * mu;
  }else{
    // the max is necessary because some combination of y and mu give negative results:
    // e.g. y = 1, mu = 0.99999999999994
    // return std::max(2.0 * (y * std::log(y/mu) - (y - mu)), 0.0);
    // However, here I don't really need this absurd precision
    return 2.0 * (y * std::log((y/mu)) - (y - mu));
  }
}


void find_closest_cluster_mb(/*Out-parameter*/ IntegerVector& cluster_assignment,
                            const NumericMatrix& Y, const NumericVector& size_factors,
                            const NumericMatrix& centers, const IntegerVector& mini_batch){
  // Rcout << "In find_closest_cluster\n";

  const int n_samples = mini_batch.length();
  const int k = centers.ncol();
  const int n_features = Y.nrow();

  for(int mbi = 0; mbi < n_samples; mbi++){
    int si = mini_batch[mbi];
    NumericMatrix::ConstColumn y = Y.column(si);
    double sf = size_factors[si];
    double min = std::numeric_limits<double>::infinity();
    int assignment = -1;
    for(int ki = 0; ki < k; ki++){
      NumericMatrix::ConstColumn mu = centers.column(ki);
      double dev = 0.0;
      for(int ri = 0; ri < n_features; ri++){
        dev += poisson_deviance_impl(y[ri], mu[ri] * sf);
      }
      if(dev < min){
        assignment = ki;
        min = dev;
      }
    }
    cluster_assignment[si] = assignment;
  }

}

void update_cluster_centers_mb(/*Out-parameter*/ NumericMatrix& centers,
                              const NumericMatrix& Y, const NumericVector& size_factors,
                              const IntegerVector& cluster_assignment,
                              const IntegerVector& mini_batch,
                              const double min_mu){
  // Rcout << "In update_cluster_centers\n";

  const int n_samples = mini_batch.length();
  const int k = centers.ncol();
  const int n_features = Y.nrow();
  centers.fill(0);
  NumericVector center_sf(k);

  for(int mbi = 0; mbi < n_samples; mbi++){
    int si = mini_batch[mbi];
    NumericMatrix::ConstColumn y = Y.column(si);
    NumericMatrix::Column col = centers.column(cluster_assignment[si]);
    center_sf[cluster_assignment[si]] += size_factors[si];
    for(int ri = 0; ri < n_features; ri++){
      col[ri] += y[ri];
    }
  }
  for(int ki = 0; ki < k; ki++){
    NumericMatrix::Column col = centers.column(ki);
    double sf = center_sf[ki];
    for(int ri = 0; ri < n_features; ri++){
      if(col[ri] == 0){
        col[ri] = min_mu;
      }else{
        col[ri] = col[ri] / sf;
      }
    }
  }
}

double calculate_error_mb(const NumericMatrix& centers, const NumericMatrix& old_centers){
  double total = 0;
  const int k = centers.ncol();
  const int n_features = centers.nrow();

  for(int ki = 0; ki < k; ki++){
    NumericMatrix::ConstColumn col = centers.column(ki);
    NumericMatrix::ConstColumn old_col = old_centers.column(ki);
    for(int ri = 0; ri < n_features; ri++){
      total += pow(col[ri] - old_col[ri], 2);
    }
  }
  return total;
}

void save_old_centers_mb(/*InOut-parameter*/NumericMatrix& a,
                      /*InOut-parameter*/NumericMatrix& b){
  NumericMatrix tmp = a;
  a = b;
  b = tmp;
}

// [[Rcpp::export]]
List run_kmeans_mb(const NumericMatrix& Y, const NumericVector& size_factors,
                   const NumericMatrix& centers_start, int mini_batch_size,
                   double min_mu, int max_iter, double tolerance,
                   bool verbose){
  // Rcout << "In run_kmeans\n";

  int n_samples = Y.ncol();
  IntegerVector cluster_assignment(n_samples);
  NumericMatrix centers = clone(centers_start);
  NumericMatrix old_centers = clone(centers);
  int iter = 0;
  for(; iter < max_iter; iter++){

    IntegerVector mini_batch = sample(n_samples, /*size=*/mini_batch_size, /*replace=*/false,
                                      /*probs=*/R_NilValue, /*one_based=*/false);

    // Changes the cluster_assignment variable
    find_closest_cluster_mb(cluster_assignment, Y, size_factors, centers, mini_batch);

    // Stores centers in old_centers and provides a matrix with correct
    // size that will be overiden anyway by update_cluster_centers
    save_old_centers_mb(centers, old_centers);

    // Changes the centers variable
    update_cluster_centers_mb(centers, Y, size_factors, cluster_assignment, mini_batch, min_mu);

    // Check for convergence
    double error = calculate_error_mb(centers, old_centers);
    if(verbose){
      Rcout << "iter " << iter << "\terror " << error << "\n";
    }
    if(error < tolerance){\
      if(verbose){
        Rcout << "Converged, after " << iter <<  " iterations. Error: " << error;
      }
      // The algorithm has converged
      break;
    }

  }


  return List::create(Named("centers") = centers,
                      Named("cluster") = cluster_assignment + 1,
                      Named("iterations") = iter + 1);

}





// // [[Rcpp::export]]
// NumericMatrix kmeans_pp_initialization(const NumericMatrix& Y, const NumericVector& size_factors,
//                                        const int k, const double min_mu){
//   const int n_samples = Y.ncol();
//   const int n_features = Y.nrow();
//
//   NumericMatrix centers(n_features, k);
//   IntegerVector center_indices(k);
//   NumericVector dev(n_samples);
//
//   int c1_idx = sample(n_samples, /*size=*/1, /*replace=*/false, /*probs=*/R_NilValue, /*one_based=*/false)[0];
//   centers.column(0) = pmax(Y.column(c1_idx) / size_factors[c1_idx], min_mu);
//   NumericMatrix::Column c1 = centers.column(0);
//   for(int si = 0; si < n_samples; si++){
//     double dev_tmp = 0.0;
//     NumericMatrix::ConstColumn y = Y.column(si);
//     for(int ri = 0; ri < n_features; ri++){
//       dev_tmp += poisson_deviance_impl(y[ri], c1[ri]);
//     }
//     dev[si] = dev_tmp;
//   }
//
//   int ci_idx = -1;
//   for(int ki = 1; ki < k; ki++){
//     // Draw random number until it is different from values in center_indices
//     bool search_again = true;
//     while(search_again){
//       // Rcout << "Search for next index\n" << dev << "\n";
//       ci_idx = sample(n_samples, /*size=*/1, /*replace=*/false, /*probs=*/dev,  /*one_based=*/false)[0];
//       search_again = false;
//       for(int kii = 0; kii < ki; kii++){
//         if(center_indices[kii] == ci_idx){
//           search_again = true;
//           break;
//         }
//       }
//     }
//
//     centers.column(ki)= pmax(Y.column(c1_idx) / size_factors[c1_idx], min_mu);
//     NumericMatrix::Column ci = centers.column(ki);
//
//     if(ki == k){
//       // Save the last dev recalculation
//       break;
//     }
//
//     for(int si = 0; si < n_samples; si++){
//       double dev_tmp = 0.0;
//       NumericMatrix::ConstColumn y = Y.column(si);
//       for(int ri = 0; ri < n_features; ri++){
//         dev_tmp += poisson_deviance_impl(y[ri], ci[ri]);
//       }
//       dev[si] = dev_tmp;
//     }
//   }
//   return centers;
// }
//
//
//
// // [[Rcpp::export(rng=false)]]
// double benchmark_log(const NumericVector&  y, const NumericVector&  mu){
//   int n = y.length();
//   double total = 0;
//   for(int i = 0; i < n; i++){
//     total += poisson_deviance_impl(y[i], mu[i]);
//   }
//   return total;
// }

