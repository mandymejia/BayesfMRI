#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
// typedef Eigen::Triplet<double> T;

//' Update the value of kappa2
//'
//' @param phi a scalar
//' @param in_list a list with elements Cmat, Gmat, and GtCinvG
//' @param n_sess the integer number of sessions
//' @param a_star pre-computed coefficient
//' @param b_star pre-computed coefficient
//' @param tol the numeric tolerance for finding the optimal value of kappa2
//' @export
// [[Rcpp::export]]
double updateKappa2(double phi, Rcpp::List in_list, int n_sess, double a_star, double b_star, double tol) {
  // Read in the SPDE matrices
  SpMat Cmat     = SpMat (in_list["Cmat"]);
  SpMat Gmat     = SpMat (in_list["Gmat"]);
  SpMat GtCinvG     = SpMat (in_list["GtCinvG"]);

  // Initialize the Q matrix
  SpMat Q = Cmat + 2 * Gmat + GtCinvG;
  int p = Q.rows();

  // Find the symbolic Cholesky factorization of Q
  SimplicialLLT<SpMat> cholQ(Q);
  cholQ.analyzePattern(Q);
  SpMat L = cholQ.matrixL();
  // Initialize objects to find the log determinant of Q
  VectorXd diagL(p);
  VectorXd lnDiagL(p);
  // Grab diagonals from the Cholesky decomposition and take their log
  for (int i = 0; i < p; ++i) {
    if(L.coeff(i,i) != 0) diagL(i) = L.coeff(i,i);
    lnDiagL(i) = log(diagL(i));
  }
  // Find the log of the determinant of Q
  double lDQ = n_sess * 2 * lnDiagL.sum();

  // Evaluate objective function on a grid from 0 to 50
  int length_grid = 50 / tol;
  // Initialize the k_out and objFn values
  double objFn_old = a_star + b_star - lDQ;
  double objFn = 0;
  double k_star = 0;
  double k_out = 0;
  // Work through the grid and find the value of kappa2 that minimizes the objFn
  for(int l = 0; l < length_grid; ++l) {
    // Increment k_star
    k_star += tol;
    // Calculate Q
    Q = k_star * Cmat + 2 * Gmat + GtCinvG / k_star;
    // Use the symbolic Cholesky to update the Cholesky of Q
    cholQ.factorize(Q);
    // Take the lower triangular Cholesky matrix
    L = cholQ.matrixL();
    // Grab diagonals from the Cholesky decomposition and take their log
    for (int i = 0; i < p; ++i) {
      if(L.coeff(i,i) != 0) diagL(i) = L.coeff(i,i);
      lnDiagL(i) = log(diagL(i));
    }
    // Sum the logged values, then multiply by 2 and the number of sessions
    lDQ = n_sess * 2 * lnDiagL.sum();
    // Find the value of the objective function
    objFn = k_star * a_star + b_star / k_star - lDQ;
    // If the objective function is less than any previous value, update k_out
    if(objFn < objFn_old) {
      k_out = k_star;
      objFn_old = objFn;
    }
  }
  return k_out;
}

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

