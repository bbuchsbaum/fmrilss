library(fmrilss)

# Create minimal test data
set.seed(123)
n <- 10
T <- 5
V <- 3

C <- matrix(rnorm(n*T), n, T)
Y <- matrix(rnorm(n*V), n, V)

cat("Testing C++ functions step by step...\n")
cat("C dims:", dim(C), "\n")
cat("Y dims:", dim(Y), "\n")

# Test if basic operations work
cat("Testing basic operations...\n")
tryCatch({
  sumC <- sum(C, 1)
  cat("sum(C, 1) - result type:", class(sumC), "dims:", length(sumC), "\n")
}, error = function(e) {
  cat("Error in sum(C, 1):", e$message, "\n")
})

# Create a simple C++ test function to debug step by step
cat("Testing C++ dimension step by step...\n")
cpp_debug_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int debug_dimensions(const arma::mat& C, const arma::mat& Y) {
  
  Rcout << "C dimensions: " << C.n_rows << "x" << C.n_cols << std::endl;
  Rcout << "Y dimensions: " << Y.n_rows << "x" << Y.n_cols << std::endl;
  
  // Test what sum(C, 1) actually returns  
  auto sumC_result = sum(C, 1);
  Rcout << "sum(C,1) type: " << typeid(sumC_result).name() << std::endl;
  Rcout << "sum(C,1) n_rows: " << sumC_result.n_rows << " n_cols: " << sumC_result.n_cols << std::endl;
  
  // Try to convert to rowvec
  try {
    arma::rowvec sumC = sumC_result;
    Rcout << "Conversion to rowvec successful" << std::endl;
    Rcout << "sumC dims after conversion: " << sumC.n_rows << "x" << sumC.n_cols << std::endl;
    return 1;
  } catch(const std::exception& e) {
    Rcout << "Error converting to rowvec: " << e.what() << std::endl;
    return 0;
  }
}
'

# Compile and run the debug function
library(Rcpp)
tryCatch({
  sourceCpp(code = cpp_debug_code)
  result <- debug_dimensions(C, Y)
  cat("Debug completed with result:", result, "\n")
}, error = function(e) {
  cat("Error in C++ debug:", e$message, "\n")
})

# Test C++ function call
cat("Testing lss_compute_cpp...\n")
tryCatch({
  result <- fmrilss:::lss_compute_cpp(C, Y)
  cat("Success! Result dims:", dim(result), "\n")
}, error = function(e) {
  cat("Error in lss_compute_cpp:", e$message, "\n")
}) 