#include <algorithm>    // std::min
#include <Rcpp.h>
using namespace Rcpp;



// The log_approx implementation relies on a lookup table for a reduced space of floating point numbers
// I use 6 bits for the exponent and 8 bits for the significand. This totals to a lookup vector with
// 2^6 * 2^8 = 2^14 or 32 kilobyte (because entries are doubles)
// The vector is called LOG_LOOKUP_VECTOR
const int LOG_LOOKUP_TABLE_EXPONENTIAL_BITS = 8;
const int LOG_LOOKUP_TABLE_SIGNIFCAND_BITS = 6;
const int LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH = pow(2, LOG_LOOKUP_TABLE_EXPONENTIAL_BITS);
const int LOG_LOOKUP_TABLE_EXPONENTIAL_OFFSET = LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH / 2;
const int LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH = pow(2, LOG_LOOKUP_TABLE_SIGNIFCAND_BITS);
const int LOG_LOOKUP_TABLE_LENGTH = LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH * LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH;
const int LOG_LOOKUP_TABLE_LENGTH_MINUS_ONE = LOG_LOOKUP_TABLE_LENGTH - 1;

const NumericVector LOG_LOOUP_VECTOR = (NumericVector) Environment::namespace_env("kmeansPoi")["LOG_LOOKUP_VECTOR"];
const NumericVector LOG_LOOUP_VECTOR2 = (NumericVector) Environment::namespace_env("kmeansPoi")["LOG_LOOKUP_VECTOR2"];

// [[Rcpp::export]]
double fuse_ints_for_double(int exp, int sig){
  long long tmp = ((long)exp + 1023) << 52;
  tmp = tmp | (((long)sig) << (52 - LOG_LOOKUP_TABLE_SIGNIFCAND_BITS));
  return *((double*) &tmp);
}

inline int sign_of_int(int x){
  return ((x < 0) * -2) + 1;
}

inline int get_exponent_as_integer_from_double_impl(double x){
  long long tmp = *(long long*) &x;
  return (int) ((tmp & 0x7FF0000000000000) >> 52) - 1023;
}

inline int get_significand_as_integer_from_double_impl(double x){
  long long tmp = *(long long*) &x;
  return (int) ((tmp & 0xFFFFFFFFFFFFF) >> (52 - LOG_LOOKUP_TABLE_SIGNIFCAND_BITS));
}

inline double log_approx_impl2(const double x){
  int index = (get_exponent_as_integer_from_double_impl(x)+LOG_LOOKUP_TABLE_EXPONENTIAL_OFFSET) * LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH + get_significand_as_integer_from_double_impl(x);
  index = index < 0 ? 0 : index >= LOG_LOOKUP_TABLE_LENGTH ? LOG_LOOKUP_TABLE_LENGTH_MINUS_ONE : index;
  return LOG_LOOUP_VECTOR2[index];
}



// [[Rcpp::export]]
int get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH(){
  return LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH;
}
// [[Rcpp::export]]
int get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH(){
  return LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH;
}


// [[Rcpp::export]]
int get_exponent_as_integer_from_double(double x){
  return get_exponent_as_integer_from_double_impl(x);
}

// [[Rcpp::export]]
int get_significand_as_integer_from_double(double x){
  return get_significand_as_integer_from_double_impl(x);
}

// [[Rcpp::export]]
NumericVector log_approx2(const NumericVector& x){
  int n = x.length();
  NumericVector result(n);
  for(int i = 0; i < n; i++){
    result[i] = log_approx_impl2(x[i]);
  }
  return result;
}




inline double poisson_deviance_impl2(const double y, const double mu){
  return 2.0 * (y * (log_approx_impl2(y/mu) - 1) + mu);
}


// [[Rcpp::export]]
NumericVector poisson_deviance2(const NumericVector& y, const NumericVector& mu){
  if(y.length() == mu.length()){
    int max_iter = y.length();
    NumericVector result(y.length());
    for(int i = 0; i < max_iter; i++){
      result[i] = poisson_deviance_impl2(y[i], mu[i]);
    }
    return result;
  }else if(y.length() == 1){
    int max_iter = mu.length();
    NumericVector result(max_iter);
    for(int i = 0; i < max_iter; i++){
      result[i] = poisson_deviance_impl2(y[0], mu[i]);
    }
    return result;
  }else if(mu.length() == 1){
    int max_iter = y.length();
    NumericVector result(max_iter);
    for(int i = 0; i < max_iter; i++){
      result[i] = poisson_deviance_impl2(y[i], mu[0]);
    }
    return result;
  }else{
    stop("Length of y and mu differ");
  }
}


double increment_deviance(const NumericVector& y, const NumericVector& mu, const int n_features, const double sf){
  double dev = 0.0;
  for(int ri = 0; ri < n_features; ri++){
    dev += poisson_deviance_impl2(y[ri], mu[ri] * sf);
  }
  return dev;
}



void find_closest_cluster2(/*Out-parameter*/ IntegerVector& cluster_assignment,
                          const NumericMatrix& Y, const NumericVector& size_factors,
                          const NumericMatrix& centers){
  // Rcout << "In find_closest_cluster\n";

  const int n_samples = Y.ncol();
  const int k = centers.ncol();
  const int n_features = Y.nrow();

  for(int si = 0; si < n_samples; si++){
    double sf = size_factors[si];
    double min = std::numeric_limits<double>::infinity();
    int assignment = -1;
    for(int ki = 0; ki < k; ki++){
      double dev = 0.0;
      for(int ri = 0; ri < n_features; ri++){
        dev += poisson_deviance_impl2(Y[si * n_features + ri], centers[ki * n_features + ri] * sf);
      }
      if(dev <= min){
        assignment = ki;
        min = dev;
      }
    }
    cluster_assignment[si] = assignment;
  }

}

void update_cluster_centers2(/*Out-parameter*/ NumericMatrix& centers,
                            const NumericMatrix& Y, const NumericVector& size_factors,
                            const IntegerVector& cluster_assignment,
                            const double min_mu){
  // Rcout << "In update_cluster_centers\n";

  const int n_samples = Y.ncol();
  const int k = centers.ncol();
  const int n_features = Y.nrow();
  centers.fill(0);
  NumericVector center_sf(k);

  for(int si = 0; si < n_samples; si++){
    int ca = cluster_assignment[si];
    center_sf[ca] += size_factors[si];
    for(int ri = 0; ri < n_features; ri++){
      centers[ca * n_features + ri] += Y[si * n_features + ri];
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

double calculate_error2(const NumericMatrix& centers, const NumericMatrix& old_centers){
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

void save_old_centers2(/*InOut-parameter*/NumericMatrix& a,
                      /*InOut-parameter*/NumericMatrix& b){
  NumericMatrix tmp = a;
  a = b;
  b = tmp;
}

// [[Rcpp::export]]
List run_kmeans2(const NumericMatrix& Y, const NumericVector& size_factors,
                const NumericMatrix& centers_start,
                double min_mu, int max_iter, double tolerance,
                bool verbose){
  // Rcout << "In run_kmeans\n";

  IntegerVector cluster_assignment(Y.ncol());
  NumericMatrix centers = clone(centers_start);
  NumericMatrix old_centers = clone(centers);
  int iter = 0;
  for(; iter < max_iter; iter++){

    // Rcout << "find_closest_cluster2\n";
    // Changes the cluster_assignment variable
    find_closest_cluster2(cluster_assignment, Y, size_factors, centers);

    // Rcout << "save_old_centers2\n";
    // Stores centers in old_centers and provides a matrix with correct
    // size that will be overiden anyway by update_cluster_centers
    save_old_centers2(centers, old_centers);

    // Rcout << "update_cluster_centers2\n";
    // Changes the centers variable
    update_cluster_centers2(centers, Y, size_factors, cluster_assignment, min_mu);

    // Rcout << "calculate_error2\n";
    // Check for convergence
    double error = calculate_error2(centers, old_centers);
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



// [[Rcpp::export]]
NumericMatrix kmeans_pp_initialization2(const NumericMatrix& Y, const NumericVector& size_factors,
                                        const int k, const double min_mu){
  const int n_samples = Y.ncol();
  const int n_features = Y.nrow();

  NumericMatrix centers(n_features, k);
  IntegerVector center_indices(k);
  NumericVector dev(n_samples);

  int c1_idx = sample(n_samples, /*size=*/1, /*replace=*/false, /*probs=*/R_NilValue, /*one_based=*/false)[0];
  for(int ri = 0; ri < n_features; ri++){
    centers[ri] = std::max(Y[c1_idx * n_features + ri] / size_factors[c1_idx], min_mu);
  }
  for(int si = 0; si < n_samples; si++){
    double dev_tmp = 0.0;
    for(int ri = 0; ri < n_features; ri++){
      dev_tmp += poisson_deviance_impl2(Y[si * n_features + ri], centers[ri]);
    }
    dev[si] = std::abs(dev_tmp);
  }

  int ci_idx = -1;
  for(int ki = 1; ki < k; ki++){
    // Draw random number until it is different from values in center_indices
    bool search_again = true;
    while(search_again){
      ci_idx = sample(n_samples, /*size=*/1, /*replace=*/false, /*probs=*/dev,  /*one_based=*/false)[0];
      search_again = false;
      for(int kii = 0; kii < ki; kii++){
        if(center_indices[kii] == ci_idx){
          search_again = true;
          break;
        }
      }
    }

    for(int ri = 0; ri < n_features; ri++){
      centers[ki * n_features + ri] = std::max(Y[ci_idx * n_features + ri] / size_factors[ci_idx], min_mu);
    }

    if(ki == k){
      // Save the last dev recalculation
      break;
    }

    for(int si = 0; si < n_samples; si++){
      double dev_tmp = 0.0;
      for(int ri = 0; ri < n_features; ri++){
        dev_tmp += poisson_deviance_impl2(Y[si * n_features + ri], centers[ki * n_features + ri]);
      }
      dev[si] = std::abs(dev_tmp);
    }
  }
  return centers;
}

// [[Rcpp::export(rng=false)]]
double benchmark_log2(const NumericVector&  y, const NumericVector&  mu){
  int n = y.length();
  double total = 0;
  for(int i = 0; i < n; i++){
    total += poisson_deviance_impl2(y[i], mu[i]);
  }
  return total;
}


