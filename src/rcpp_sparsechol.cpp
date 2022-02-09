#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
// typedef Eigen::Triplet<double> T;

//' Find the log of the determinant of Q
//'
//' @param kappa2 a scalar
//' @param in_list a list with elements Cmat, Gmat, and GtCinvG
//' @param n_sess the integer number of sessions
//' @export
// [[Rcpp::export]]
double logDetQ(double kappa2, Rcpp::List in_list, int n_sess) {
  //**********************************
  //      basic parameter
  //**********************************
  // double kappa2       = Rcpp::as< double > (in_list["kappa2"]);
  // Rcpp::Rcout << "kappa2 = " << kappa2 << std::endl;
  SpMat Cmat     = SpMat (in_list["Cmat"]);
  // Rcpp::Rcout << "Cmat dims: " << Cmat.rows() << " x " << Cmat.cols() << std::endl;
  SpMat Gmat     = SpMat (in_list["Gmat"]);
  // Rcpp::Rcout << "Gmat dims: " << Gmat.rows() << " x " << Gmat.cols() << std::endl;
  SpMat GtCinvG     = SpMat (in_list["GtCinvG"]);
  // Rcpp::Rcout << "GtCinvG dims: " << GtCinvG.rows() << " x " << GtCinvG.cols() << std::endl;
  // int n_sess = in_list["n_sess"];
  // Rcpp::Rcout << "n_sess = " << n_sess << std::endl;

  //**********************************
  // Initialize
  //**********************************
  // Create SparseMatrix Q
  SpMat Qstart = kappa2 * Cmat + 2 * Gmat + GtCinvG / kappa2;
  // Rcpp::Rcout << "Qstart dims: " << Qstart.rows() << " x " << Qstart.cols() << ", num nonzero: " << Qstart.nonZeros() << std::endl;
  int p = Qstart.rows();
  // Take the Cholesky decomposition of Qstart
  SimplicialLLT<SpMat> Rstart(Qstart);
  // Rcpp::Rcout << "Rstart dims: " << Rstart.rows() << " x " << Rstart.cols() << std::endl;
  // Get the symbolic decomposition of Rstart using Qstart
  // Rstart.analyzePattern(Qstart);
  // Rstart.factorize(Qstart);
  SpMat L_out = Rstart.matrixL();

  // double L11 = L_out.coeff(1,1);
  // Rcpp::Rcout << "L[1,1] is: " << L11 << std::endl;
  // Grab diagonals from the Cholesky decomposition and take their log
  VectorXd diagR(p);
  VectorXd lnDiagR(p);
  for (int i = 0; i < p; ++i) {
    if(L_out.coeff(i,i) != 0) diagR(i) = L_out.coeff(i,i);
    lnDiagR(i) = log(diagR(i));
  }
  double lDQ = n_sess * 2 * lnDiagR.sum();
  // Rcpp::Rcout << "The log determinant of Qstart is " << lDQ << std::endl;

  // std::vector<T> tripletList;
  //
  // for (int i = 0; i < n_sess; ++i) {
  //   for (int ii = i * p; ii < i * p + p; ++ii) {
  //     for (int jj = i * p; jj < i * p + p;  ++jj) {
  //       // direct assignment
  //       // Q.insert(ii, jj) = Qstart(ii - i, jj - i);
  //
  //       // triplets ??
  //       tripletList.push_back(T(ii, jj, Qstart(ii - i, jj - i)));
  //     }
  //   }
  // }
  // SparseMatrix<double> Q(p * n_sess, p * n_sess);
  // Q.setFromTriplets(tripletList.begin(), tripletList.end);
  // Rcpp::Rcout << "Q num nonzero: " << Q.nonZeros() << std::endl;

  //**********************************
  // Begin optimization
  //**********************************
  // Create SparseMatrix Q
  // SparseMatrix<double> Q = kappa2 ...;
  // Symbolic Cholesky factorization R (with analyzePattern)
  // SimplicialLLT<SparseMatrix<double>> R;
  // R.analyzePattern(Q);
  // Update R using factorize
  // R.factorize(Q);

  // General tip: test line-by-line. Compile regularly and check and test each line's output
  // Try using Rcpp::Rcout so that the C code is outputted on the screen


  // return kappa2;
  return lDQ;
}

