#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> fast_multiply(Eigen::SparseMatrix<double> M1, Eigen::SparseMatrix<double> M2){
	Eigen::SparseMatrix<double> crosspd = M1 * M2.transpose();
	return crosspd;
}

