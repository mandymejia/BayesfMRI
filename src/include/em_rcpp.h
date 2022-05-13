#include <Rcpp.h>
#include <RcppEigen.h>

// Find the log determinant of the sparse SPDE precision Q-tilde
double logDetQt(double kappa2, const Rcpp::List &in_list, double n_sess);

// Make the sparse SPDE precision matrix Q-tilde
void makeQt(Eigen::SparseMatrix<double>* Q, double kappa2, const Rcpp::List &spde);

// The objective function for kappa^2 used in the EM
double kappa2Obj(double kappa2, const Rcpp::List &spde, double a_star, double b_star, double n_sess);

// Optimization function used to find kappa2 in the EM
double kappa2Brent(double lower, double upper, const Rcpp::List &spde, double a_star, double b_star, double n_sess);

// Make the random matrix V for the Hutchinson trace approximation
Eigen::MatrixXd makeV(int n_spde, int Ns);

// Replace sparse block matrices within existing (larger) sparse matrix
void setSparseBlock_update(Eigen::SparseMatrix<double,0,int>* A,int i, int j, Eigen::SparseMatrix<double,0,int>& B);

// Create a sparse block diagonal matrix using a list of sparse matrices
Eigen::SparseMatrix<double> sparseBdiag(Rcpp::List B_list);

//Global Control Variable
struct SquaremControl{
  int K=1;
  int method=3;//1,2,3 indicates the types of step length to be used in squarem1,squarem2, 4,5 for "rre" and "mpe" in cyclem1 and cyclem2,  standing for reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation.
  // K=1 must go with method=1,2 or 3
  // K>1 must go with method=4 or 5
  double mstep=4;
  int maxiter=1500;
  bool square=true;
  bool trace=true;//currently set to be true for debugging purpose
  double stepmin0=1;
  double stepmax0=1;
  double kr=1;
  double objfninc=1;//0 to enforce monotonicity, Inf for non-monotonic scheme, 1 for monotonicity far from solution and allows for non-monotonicity closer to solution
  double tol=1e-7;
} SquaremDefault;

//Output Struct
struct SquaremOutput{
  Eigen::VectorXd par;
  double valueobjfn;
  int iter=0;
  int pfevals=0;
  int objfevals=0;
  bool convergence=false;
} sqobj,sqobjnull;
