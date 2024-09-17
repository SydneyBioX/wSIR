// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export(.matMultArma)]]
arma::mat matMultArma(const arma::mat& A, const arma::mat& B) {
  return A * B;
}

// [[Rcpp::export(.computeRandZ)]]
Rcpp::List computeRandZ(const arma::mat& X) {
  int n = X.n_rows;

  // Center the matrix X
  arma::mat Xc = X.each_row() - arma::mean(X, 0);

  // Perform QR decomposition
  arma::mat Q, R;
  qr_econ(Q, R, Xc);

  // Adjust R to be upper triangular and ensure correct size
  R = R.submat(0, 0, Xc.n_cols - 1, Xc.n_cols - 1);

  // Multiply Q by sqrt(n)
  arma::mat Z = Q * std::sqrt(n);

  // Return R and Z as a list
  return Rcpp::List::create(Rcpp::Named("R") = R,
                            Rcpp::Named("Z") = Z);
}


// [[Rcpp::export(.fastEigen)]]
Rcpp::List fastEigen(const arma::mat& S){

  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, S);
  return Rcpp::List::create(Rcpp::Named("values") = eigval,
                            Rcpp::Named("vectors") = eigvec);
}


// [[Rcpp::export(.subsetLowerTri)]]
arma::vec subsetLowerTri(const arma::mat& m) {
  // Get the number of rows and columns
  int n_rows = m.n_rows;
  int n_cols = m.n_cols;

  // Create a vector to store the elements of the lower triangle
  std::vector<double> lower_tri_elements;

  // Loop through the lower triangular part (excluding the diagonal)
  for (int i = 1; i < n_rows; ++i) {
    for (int j = 0; j < i && j < n_cols; ++j) {
      lower_tri_elements.push_back(m(i, j));
    }
  }

  // Convert the vector to an Armadillo vector
  arma::vec result = arma::conv_to<arma::vec>::from(lower_tri_elements);

  return result;
}

// Function to compute the rank of elements in a vector

// [[Rcpp::export(.rank)]]
arma::vec rank(const arma::vec &x) {
  arma::uvec indices = arma::sort_index(x);
  arma::vec ranks(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; i++) {
    ranks(indices(i)) = i + 1;
  }
  return ranks;
}

// Function to compute the Spearman correlation

// [[Rcpp::export(.spearman_correlation)]]
arma::vec spearman_correlation(const arma::vec &x, const arma::vec &y) {
  // Ensure the vectors have the same length
  if (x.n_elem != y.n_elem) {
    Rcpp::stop("Vectors must have the same length.");
  }

  // Rank the vectors
  arma::vec rx = rank(x);
  arma::vec ry = rank(y);

  // Compute the Pearson correlation on the ranks
  arma::vec spearman_corr = arma::cor(rx, ry);
  return spearman_corr;
}
