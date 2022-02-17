#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// typedef Eigen::SparseMatrix<double> Eigen::SparseMatrix<double>; // declares a column-major sparse matrix type of double

// // A functor?
// class kappaObj
// {
// private:
//   double kappa2;
// public:
//   kappaObj(double k) : kappa2(k) { }
//
//   // This operator enables calling the
//   // operator function () on objects of kappaObj
//   double operator () (List in_list, SimplicialLLT<Eigen::SparseMatrix<double>> cholQ, double a_star, double b_star, int n_sess) const {
//     // Read in the SPDE matrices
//     Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (in_list["Cmat"]);
//     Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (in_list["Gmat"]);
//     Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (in_list["GtCinvG"]);
//     // Construct the Q matrix
//     Eigen::SparseMatrix<double> Q = kappa2*Cmat + 2 * Gmat + GtCinvG / kappa2;
//     int p = Q.rows();
//     // Use the symbolic Cholesky to update the Cholesky of Q
//     cholQ.factorize(Q);
//     // Take the lower triangular Cholesky matrix
//     Eigen::SparseMatrix<double> L = cholQ.matrixL();
//     // Grab diagonals from the Cholesky decomposition and take their log
//     double lDQ = 0;
//     for (int i = 0; i < p; ++i) {
//       lDQ += log(L.diagonal()[i]);
//     }
//     // Sum the logged values, then multiply by 2 and the number of sessions
//     lDQ = n_sess * 2 * lDQ;
//     // Find the value of the objective function
//     double objFn = kappa2 * a_star + b_star / kappa2 - lDQ;
//     return objFn;
//   }
// };

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
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (in_list["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (in_list["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (in_list["GtCinvG"]);

  // Initialize the Q matrix
  Eigen::SparseMatrix<double> Q = Cmat + 2 * Gmat + GtCinvG;
  int p = Q.rows();

  // Find the symbolic Cholesky factorization of Q
  SimplicialLLT<Eigen::SparseMatrix<double>> cholQ(Q);
  // cholQ.analyzePattern(Q);
  Eigen::SparseMatrix<double> L = cholQ.matrixL();
  // Initialize objects to find the log determinant of Q
  Eigen::VectorXd diagL = L.diagonal();
  // Grab diagonals from the Cholesky decomposition and take their log
  double lDQ = 0;
  for (int i = 0; i < p; ++i) {
    lDQ += log(diagL[i]);
  }
  // Find the log of the determinant of Q
  lDQ = n_sess * 2 * lDQ;

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
    lDQ = 0;
    for (int i = 0; i < p; ++i) {
      lDQ += log(L.diagonal()[i]);
    }
    // Sum the logged values, then multiply by 2 and the number of sessions
    lDQ = n_sess * 2 * lDQ;
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
//' @param phi a scalar
//' @param spde_list a list with elements Cmat, Gmat, and GtCinvG
//' @param n_sess the integer number of sessions
//' @export
// [[Rcpp::export]]
double logDetQ(double kappa2, double phi, List spde_list, int n_sess) {
  // Load parameters
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde_list["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde_list["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde_list["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q = kappa2 * Cmat + 2 * Gmat + GtCinvG / kappa2;
  Q = Q / (4.0 * M_PI * phi);
  // Find the Cholesky, then extract the diagonals, take their log, sum them,
  // and multiply by n_sess * 2 to get the log determinant of Q
  // SimplicialLLT<Eigen::SparseMatrix<double>> cholQ(Q);
  // Eigen::SparseMatrix<double> lowerTri = cholQ.matrixL();
  // Eigen::VectorXd cholDiags = lowerTri.diagonal();
  // double lDQ = 0.0;
  // int p = Q.rows();
  // for (int i = 0; i < p; ++i) {
  //   lDQ += log(cholDiags[i]);
  // }
  // lDQ = 2.0 * (double)n_sess * lDQ;
  // This next code was suggested by David, but seems to return inconsistent results
  SimplicialLDLT<Eigen::SparseMatrix<double>> cholQ(Q);
  double lDQ = n_sess * cholQ.vectorD().array().log().sum();
  return lDQ;
}

//' Find the log of the determinant of Q_tilde
//'
//' @param kappa2 a scalar
//' @param in_list a list with elements Cmat, Gmat, and GtCinvG
//' @param n_sess the integer number of sessions
//' @export
// [[Rcpp::export]]
double logDetQt(double kappa2, Rcpp::List in_list, int n_sess) {
  // Load parameters
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (in_list["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (in_list["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (in_list["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q= kappa2 * Cmat + 2.0 * Gmat + GtCinvG / kappa2;
  // SimplicialLLT<Eigen::SparseMatrix<double>> cholQ(Q);
  SimplicialLDLT<Eigen::SparseMatrix<double>> cholQ(Q);
  // Eigen::VectorXd diagCholQ = cholQ.vectorD();
  // Eigen::VectorXd logDiagCholQ = diagCholQ.log();
  // double sumLogDiagCholQ = logDiagCholQ.sum();
  double lDQ = n_sess * cholQ.vectorD().array().log().sum();
  // double lDQ = (double)n_sess * sumLogDiagCholQ; // This returns varying values
  // int p = Q.rows();
  // // Take the Cholesky decomposition of Q
  // Eigen::SparseMatrix<double> L_out = cholQ.matrixL();
  // // Grab diagonals from the Cholesky decomposition and take the sum of their log
  // Eigen::VectorXd DiagL = L_out.diagonal();
  // double lDQ = 0;
  // for (int i = 0; i < p; ++i) {
  //   lDQ += log(DiagL[i]);
  // }
  // Multiply by 2 because Cholesky and then by the number of sessions
  // lDQ = lDQ * 2.0 * (double)n_sess;

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

  return lDQ;
}

// [[Rcpp::export]]
double logDetQtLDLT(double kappa2, Rcpp::List in_list, int n_sess) {
  // Load parameters
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (in_list["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (in_list["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (in_list["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q= kappa2 * Cmat + 2.0 * Gmat + GtCinvG / kappa2;
  // Find the LDLT decomposition
  SimplicialLDLT<Eigen::SparseMatrix<double>> cholQ(Q);
  double lDQ = cholQ.vectorD().array().log().sum();
  lDQ = n_sess * lDQ;
  return lDQ;
}

// Attempt at writing a functor based on
// https://stackoverflow.com/questions/57315008/eigen-library-levenberg-marquardt-how-to-pass-extra-information-to-minimize
// and
// https://stackoverflow.com/questions/18509228/how-to-use-the-eigen-unsupported-levenberg-marquardt-implementation
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};

// struct my_functor : Functor<double> 	{
//   double m_b;
//   int m_num_periods;
//   vector<double> m_spot_rates;
//   my_functor(double b, int num_periods, vector<double> spot_rates) : Functor<double>(14,14),m_b(b),m_num_periods(num_periods),m_spot_rates(spot_rates) {}



// // [[Rcpp::export]]
// int initKP(Eigen::VectorXd theta, List in_list, Eigen::VectorXd w, int n_sess) {
//   Eigen::VectorXd thetaExp = theta.array().exp();
//   Rcout << "theta = " << thetaExp.transpose() << std::endl;
//
//   kp_functor functor(in_list, w, n_sess);
//   // Eigen::LevenbergMarquardt<kp_functor, double> lm(kp_functor);
//   Eigen::NumericalDiff<kp_functor> numDiff(functor);
//   Eigen::LevenbergMarquardt<Eigen::NumericalDiff<kp_functor>,double> lm(numDiff);
//   lm.parameters.maxfev = 2000;
//   lm.parameters.xtol = 1.0e-5;
//
//   lm.minimize(theta);
//   // int ret = lm.minimize(theta);
//   Rcout << "lm.fjac: " << lm.fjac << std::endl;
//   Rcout << "Number of iterations: " << lm.iter << std::endl;
//   // Rcout << "Whatever ret is: " << ret << std::endl;
//   // thetaExp = exp(theta);
//   std::cout << "theta that minimizes the function: " << theta.array().exp() << std::endl;
//   return 0;
// }

struct kinit_functor : Functor<double>
{
  double m_phi;
  List m_spde;
  Eigen::VectorXd m_w;
  int m_n_sess;
  kinit_functor(double phi, List spde, Eigen::VectorXd w, int n_sess) : Functor<double>(1,1),m_phi(phi),m_spde(spde),m_w(w),m_n_sess(n_sess) {}
  int operator()(const Eigen::VectorXd &kappa2, Eigen::VectorXd &fvec) const
  {
    Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (m_spde["Cmat"]);
    Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (m_spde["Gmat"]);
    Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (m_spde["GtCinvG"]);
    // Create SparseMatrix Q
    Eigen::SparseMatrix<double> Q = kappa2(0) * Cmat + 2.0 * Gmat + GtCinvG / kappa2(0);
    Q = Q / (4.0 * M_PI);
    double lDQ = logDetQt(kappa2(0), m_spde, m_n_sess);
    Eigen::VectorXd Qw =  Q * m_w;
    double wQw = m_w.transpose() * Qw;
    fvec(0) = wQw / (8.0 * M_PI * m_phi) - lDQ;
    // fvec(1) = -lDQ;
    return 0;
  }

  // int df(const Eigen::VectorXd &kappa2, Eigen::MatrixXd &fjac) const
  // {
  //   Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (m_spde["Cmat"]);
  //   Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (m_spde["Gmat"]);
  //   Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (m_spde["GtCinvG"]);
  //   // Create SparseMatrix Q_tilde
  //   Eigen::SparseMatrix<double> Qt = kappa2(0) * Cmat + 2.0 * Gmat + GtCinvG / kappa2(0);
  //   SimplicialLLT<Eigen::SparseMatrix<double>> cholQt(Qt);
  //   Eigen::SparseMatrix<double> dQdK = Cmat - GtCinvG / (kappa2(0) * kappa2(0));
  //   Eigen::SparseMatrix<double> QtInv = cholQt.solve(Qt);
  //   Eigen::SparseMatrix<double> QtInv_dQdK = QtInv.cwiseProduct(dQdK);
  //   double first_term = QtInv_dQdK.sum();
  //   double wGCGw = m_w.transpose() * GtCinvG * m_w;
  //   double wCw = m_w.transpose() * Cmat * m_w;
  //   double second_term = wGCGw / (kappa2(0) * kappa2(0)) - wCw;
  //   second_term = second_term / (4.0 * M_PI * m_phi);
  //
  //   fjac(0,0) = first_term + second_term;
  //
  //   return 0;
  // }
};

// [[Rcpp::export]]
Eigen::VectorXd initK(Eigen::VectorXd kappa2, double phi, List in_list, Eigen::VectorXd w, int n_sess) {
  // Rcout << "kappa2 init = " << kappa2 << std::endl;

  kinit_functor functor(phi, in_list, w, n_sess);
  // LevenbergMarquardt<k_functor, double> lm(functor);
  Eigen::NumericalDiff<kinit_functor> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<kinit_functor>,double> lm(numDiff);
  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1.0e-5;

  lm.minimize(kappa2);
  // int ret = lm.minimize(kappa2);
  // Rcout << "lm.fjac: " << lm.fjac << std::endl;
  // Rcout << "Number of iterations: " << lm.iter << std::endl;
  // Rcout << "Whatever ret is: " << ret << std::endl;

  // std::cout << "kappa2 that minimizes the function: " << kappa2 << std::endl
  return kappa2;
}

struct updateKfunctor : Functor<double>
{
  List m_spde;
  double m_a_star;
  double m_b_star;
  int m_n_sess;
  updateKfunctor(List spde, double a_star, double b_star, int n_sess) : Functor<double>(1,1),m_spde(spde),m_a_star(a_star),m_b_star(b_star),m_n_sess(n_sess) {}
  int operator()(const VectorXd &kappa2, Eigen::VectorXd &fvec) const
  {
    VectorXd expKappa2 = kappa2.exp();
    double lDQ = logDetQt(expKappa2(0), m_spde, m_n_sess);
    fvec(0) = expKappa2(0) * m_a_star + m_b_star / expKappa2(0) - lDQ;
    // fvec(1) = -lDQ;
    return 0;
  }
};

// [[Rcpp::export]]
Eigen::VectorXd updateK(VectorXd kappa2, List spde, double a_star, double b_star, int n_sess) {
  updateKfunctor functor(spde, a_star, b_star, n_sess);
  Eigen::NumericalDiff<updateKfunctor> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<updateKfunctor>,double> lm(numDiff);
  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1.0e-5;
  lm.minimize(kappa2);
  return kappa2;
}

//' Find the initial values of kappa2 and phi
//'
//' @param kappa2 a scalar scale parameter
//' @param phi a scalar range parameter
//' @param spde a list containing the sparse matrix elements Cmat, Gmat, and GtCinvG
//' @param w the beta_hat estimates for a single task
//' @param n_sess the number of sessions
//' @param tol the stopping rule tolerance
//' @export
// [[Rcpp::export]]
Eigen::VectorXd initialKP(Eigen::VectorXd kappa2, double phi, List spde, Eigen::VectorXd w, int n_sess, double tol) {
  Eigen::VectorXd old_theta(2);
  old_theta(0) = kappa2(0);
  old_theta(1) = phi;
  kappa2 = initK(kappa2, phi, spde, w, n_sess);
  // Rcout << "First kappa2 is " << kappa2 << std::endl;
  int n_spde = w.size();
  // Rcout << "n_spde = " << n_spde << std::endl;
  n_spde = n_spde / n_sess;
  // Rcout << "n_spde / n_sess = " << n_spde << std::endl;
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q = kappa2(0) * Cmat + 2.0 * Gmat + GtCinvG / kappa2(0);
  Eigen::VectorXd Qw = Q * w;
  double wQw = w.transpose() * Qw;
  phi = wQw / (4.0 * M_PI * n_spde * n_sess);
  // Rcout << "First phi = " << phi << std::endl;
  Eigen::VectorXd theta(2);
  theta(0) = kappa2(0);
  theta(1) = phi;
  Eigen::VectorXd my_diff = theta - old_theta;
  double eps = my_diff.squaredNorm();
  // Rcout << "First epsilon = " << eps << std::endl;
  int num_steps = 1;

  while(eps > tol) {
    old_theta(0) = kappa2(0);
    old_theta(1) = phi;
    kappa2 = initK(kappa2, phi, spde, w, n_sess);
    Q = kappa2(0) * Cmat + 2.0 * Gmat + GtCinvG / kappa2(0);
    Qw = Q * w;
    wQw = w.transpose() * Qw;
    phi = wQw / (4.0 * M_PI * n_spde * n_sess);
    theta(0) = kappa2(0);
    theta(1) = phi;
    my_diff = theta - old_theta;
    eps = my_diff.squaredNorm();
    num_steps += 1;
  }
  // Rcout << "Number of steps to converge: " << num_steps << std::endl;
  return theta;
}

/*
 B of size n1 x n2
 Set A(i:i+n1,j:j+n2) = B
 */
void setSparseBlock(SparseMatrix<double,0,int>* A,int i, int j, SparseMatrix<double,0,int>& B)
{
  for (int k=0; k<B.outerSize(); ++k) {
    for (SparseMatrix<double,0,int>::InnerIterator it(B,k); it; ++it) {
      A->insert(it.row()+i, it.col()+j) = it.value();
    }
  }
}

/*
 B of size n1 x n2
 Set A(i:i+n1,j:j+n2) = B (update)
 */
void setSparseBlock_update(SparseMatrix<double,0,int>* A,int i, int j, SparseMatrix<double,0,int>& B)
{
  for (int k=0; k<B.outerSize(); ++k) {
    for (SparseMatrix<double,0,int>::InnerIterator it(B,k); it; ++it) {
      A->coeffRef(it.row()+i, it.col()+j) = it.value();
    }
  }
}

// [[Rcpp::export]]
Eigen::VectorXd findTheta(Eigen::VectorXd theta, List spde, Eigen::VectorXd y, Eigen::SparseMatrix<double> X, Eigen::SparseMatrix<double> QK, Eigen::SparseMatrix<double> Psi, Eigen::SparseMatrix<double> A, Eigen::MatrixXd Vh, double tol) {
  // Bring in the spde matrices
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  int Ns = 50; // For the Hutchinson estimator involving the trace of matrix products with Sigma
  int n_spde = Cmat.rows();
  int nKs = A.rows();
  int K = theta.size();
  K = (K - 1) / 2;
  int n_sess = nKs / (n_spde * K);
  int sig2_ind = 2*K;
  SimplicialLLT<Eigen::SparseMatrix<double>> cholQK(QK);
  cholQK.analyzePattern(QK);
  Eigen::VectorXd theta_old = theta;
  int Step = 0;
  Rcout << "Step 0 theta: " << theta_old.transpose() << std::endl;
  double eps = tol + 1;
  // Create SparseMatrix Q
  // SparseMatrix<double> Q = kappa2 ...;
  // Symbolic Cholesky factorization R (with analyzePattern)
  // SimplicialLLT<SparseMatrix<double>> R;
  // R.analyzePattern(Q);
  // Update R using factorize
  // R.factorize(Q);
  while(eps > tol) {
    cholQK.factorize(QK);
    Eigen::SparseMatrix<double> AdivS2 = A / theta[sig2_ind];
    Eigen::SparseMatrix<double> Sig_inv = QK + AdivS2;
    Eigen::SparseMatrix<double> Xpsi = X * Psi;
    Eigen::VectorXd XpsiY = Xpsi.transpose() * y;
    Eigen::VectorXd m = XpsiY / theta(sig2_ind);
    Eigen::VectorXd mu = cholQK.solve(m);
    // Solve for sigma_2
    Eigen::VectorXd XpsiMu = Xpsi * mu;
    Eigen::MatrixXd P = cholQK.solve(Vh);
    Eigen::MatrixXd Avh = A * Vh;
    Eigen::MatrixXd PaVh = P.transpose() * Avh;
    Eigen::VectorXd diagPaVh = PaVh.diagonal();
    double TrSigA = diagPaVh.sum() / Ns;
    Eigen::VectorXd Amu = A * mu;
    double muAmu = mu.transpose() * Amu;
    double TrAEww = muAmu + TrSigA;
    double yy = y.transpose() * y;
    double yXpsiMu = y.transpose() * XpsiMu;
    int ySize = y.size();
    theta[sig2_ind] = (yy - 2 * yXpsiMu + TrAEww) / ySize;
    // Update kappa2 and phi by task
    double a_star;
    double b_star;
    int idx_start;
    int idx_stop;
    Eigen::VectorXd muKns(n_spde);
    Eigen::VectorXd Cmu(n_spde);
    double muCmu;
    Eigen::VectorXd Gmu(n_spde);
    double muGmu;
    Eigen::MatrixXd Pkn(n_spde, Ns);
    Eigen::MatrixXd Vkn(n_spde, Ns);
    Eigen::MatrixXd CVkn(n_spde, Ns);
    Eigen::MatrixXd PCVkn(Ns, Ns);
    Eigen::VectorXd diagPCVkn(Ns);
    double sumDiagPCVkn;
    Eigen::MatrixXd GVkn(n_spde, Ns);
    Eigen::MatrixXd PGVkn(Ns, Ns);
    Eigen::VectorXd diagPGVkn(Ns);
    double sumDiagPGVkn;
    Eigen::VectorXd GCGmu(n_spde);
    double muGCGmu;
    Eigen::MatrixXd GCGVkn(n_spde, Ns);
    Eigen::MatrixXd PGCGVkn(Ns,Ns);
    Eigen::VectorXd diagPGCGVkn(Ns);
    double sumDiagPGCGVkn;
    double phi_partA;
    double phi_partB;
    double phi_partC;
    for(int k = 0; k < K; k++) {
      a_star = 0.0;
      b_star = 0.0;
      muCmu = 0.0;
      muGmu = 0.0;
      muGCGmu = 0.0;
      sumDiagPCVkn = 0.0;
      sumDiagPGVkn = 0.0;
      sumDiagPGCGVkn = 0.0;
      for(int ns = 0; ns < n_sess; ns++) {
        idx_start = k * n_spde + ns * K * n_spde;
        idx_stop = idx_start + n_spde;
        muKns = mu.segment(idx_start,idx_stop);
        Cmu = Cmat * muKns;
        muCmu += muKns.transpose() * Cmu;
        Gmu = Gmat * muKns;
        muGmu += muKns.transpose() * Gmu;
        Pkn = P.block(idx_start, 0, n_spde, Ns);
        Vkn = Vh.block(idx_start, 0, n_spde, Ns);
        CVkn = Cmat * Vkn;
        PCVkn = Pkn.transpose() * CVkn;
        diagPCVkn = PCVkn.diagonal();
        sumDiagPCVkn += diagPCVkn.sum();
        GVkn = Gmat * Vkn;
        PGVkn = Pkn.transpose() * GVkn;
        diagPGVkn = PGVkn.diagonal();
        sumDiagPGVkn += diagPGVkn.sum();
        GCGmu = GtCinvG * muKns;
        muGCGmu += muKns.transpose() * GCGmu;
        GCGVkn = GtCinvG * Vkn;
        PGCGVkn = Pkn.transpose() * GCGVkn;
        diagPGCGVkn = PGCGVkn.diagonal();
        sumDiagPGCGVkn += diagPGCGVkn.sum();
      }
      sumDiagPCVkn = sumDiagPCVkn / Ns;
      // Update kappa2
      a_star = (muCmu + sumDiagPCVkn) / (4.0 * M_PI * theta[k + K]);
      sumDiagPGCGVkn = sumDiagPGCGVkn / Ns;
      b_star = (muGCGmu + sumDiagPGCGVkn) / (4.0 * M_PI * theta[k + K]);
      VectorXd new_kappa2 = theta.segment(k,k+1);
      new_kappa2 = new_kappa2.log();
      new_kappa2 = updateK(new_kappa2, spde, a_star, b_star, n_sess);
      new_kappa2 = new_kappa2.exp();
      theta[k] = new_kappa2(0);
      // Update phi
      phi_partA = sumDiagPCVkn + muCmu;
      phi_partA = phi_partA * new_kappa2(0);
      Rcout << "phi_partA = " << phi_partA;
      phi_partB = sumDiagPGVkn + muGmu;
      phi_partB = 2 * phi_partB;
      Rcout << ", phi_partB = " << phi_partB;
      phi_partC = sumDiagPGCGVkn + muGCGmu;
      phi_partC = phi_partC / new_kappa2(0);
      Rcout << ", phi_partC = " << phi_partC << std::endl;
      double phi_new = phi_partA + phi_partB + phi_partC;
      phi_new = phi_new / (4 * M_PI * n_spde * n_sess);
      theta[k + K] = phi_new;
      Eigen::SparseMatrix<double> Qk = theta(k) * Cmat + 2.0 * Gmat + GtCinvG / theta(k);
      Qk = Qk / (4 * M_PI * theta(k + K));
      for(int ns = 0; ns < n_sess; ns++) {
        int start_i = k * n_spde + ns * K * n_spde;
        setSparseBlock_update(&QK, start_i, start_i, Qk);
      }
    }
    Step += 1;
    Rcout << "Step " << Step << " theta: " << theta.transpose() << std::endl;
    VectorXd thetaDiff = theta - theta_old;
    eps = thetaDiff.squaredNorm();
    // Rcout << "The squared difference is " << eps << std::endl;
    theta_old = theta;
  }
  return theta;
}

// // [[Rcpp::export]]
// double initK(double phi, List spde, Eigen::VectorXd w, int n_sess, double tol) {
//     Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
//     Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
//     Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
//     // Create SparseMatrix Q_tilde
//     Eigen::SparseMatrix<double> Qt = Cmat + 2.0 * Gmat + GtCinvG;
//     SimplicialLDLT<Eigen::SparseMatrix<double>> ldltQt(Qt);
//     ldltQt.analyzePattern(Qt);
//     // Set initial search triplet for the golden section search
//     Eigen::VectorXd k_sec(3);
//     k_sec(0) = 1e-5;
//     k_sec(1) = 50.0 / 1.618;
//     k_sec(2) = 50.0;
//     Eigen::VectorXd f_sec(3);
//     // double Xa = k_sec[1] - k_sec[0];
//     // double Xb = k_sec[2] - k_sec[1];
//     // Evaluate at the initial points
//     for(int i = 0; i < 3; i++){
//       Qt = k_sec[i] * Cmat + 2.0 * Gmat + GtCinvG / k_sec[i];
//       ldltQt.factorize(Qt);
//       double lDQ = n_sess * ldltQt.vectorD().array().log().sum();
//       Rcout << "i = " << i << " and lDQ = " <<  lDQ;
//       Eigen::VectorXd Qw = Qt * w;
//       double wQw = w.transpose() * Qw;
//       // wQw = wQw / (4 * M_PI* phi);
//       Rcout << " and wQw = " << wQw << std::endl;
//       f_sec[i] = wQw - lDQ;
//     }
//     Rcout << "The initial points are " << k_sec.transpose() << std::endl;
//     Rcout << "The log likelihoods are " << f_sec.transpose() << std::endl;
//   return 0.0;
// }

// struct my_functor : Functor<double>
// {
//   my_functor(void) : Functor<double>(14,14) {}
//   int operator()(const Eigen::Eigen::VectorXd &x, Eigen::Eigen::VectorXd &fvec) const
//   {
//     // These are the extra parameters I need to pass
//     //which are used in generating fvec values
//     double b = 0.005;
//     int num_periods = 14;
//     vector<double> spot_rates = { 0.073,0.0762,0.081,0.0845,0.092,0.0964,0.1012,0.1045,0.1075,0.1122,0.1155,0.1192,0.122,0.1232 };
//     //End of extra parameters
//
//     vector<double> vec(x.data(), x.data() + x.rows() * x.cols());
//     vector<vector<double>> short_rates = bdt_short_rate_lattice(spot_rates[0], vec, b, num_periods);
//     vector<vector<double>> bdt_elementary_prices = bdt_elementary_lattice(short_rates, 0.5);
//     vector<double> bdt_zcb_prices;
//     vector<double> bdt_spot_rates;
//     for (int i = 1; i <= short_rates.size(); ++i) {
//       bdt_zcb_prices.push_back(accumulate(bdt_elementary_prices[i].begin(), bdt_elementary_prices[i].end(), 0.0));
//       bdt_spot_rates.push_back(pow((1.0 / bdt_zcb_prices.back()), (1.0 / i)) - 1.0);
//     }
//     for (int i = 0; i < bdt_spot_rates.size(); ++i) {
//       fvec(i) = spot_rates[i] - bdt_spot_rates[i];
//     }
//
//     return 0;
//   }
// };


// // functor for the initial conditions problem
// struct initFunctor
// {
//   // Compute "error" that we try to minimize
//   int operator()(const VectorXf &theta, VectorXf &fvec, List spde_list, int n_sess, VectorXf w) const
//   {
//     float kappa2 = theta(0);
//     float phi = theta(1);
//     Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde_list["Cmat"]);
//     Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde_list["Gmat"]);
//     Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde_list["GtCinvG"]);
//     // Create SparseMatrix Q
//     Eigen::SparseMatrix<double> Q = kappa2 * Cmat + 2 * Gmat + GtCinvG / kappa2;
//     Q = Q / (4 * M_PI * phi);
//     float lDQ = logDetQ(kappa2, phi, spde_list, n_sess);
//     float wQw = w.transpose() * Q * w;
//     fvec(0) = wQw - lDQ;
//     fvec(1) = 0;
//     return 0;
//   }
//
//   // Compute the Jacobian of the function
//   int df(const VectorXf &theta, Eigen::MatrixXd &fjac, List spde_list, int n_sess, VectorXf w) const
//   {
//     float epsilon;
//     epsilon = 1e-5f;
//
//     for (int i = 0; i < theta.size(); ++i) {
//       // Add and subtract just a little bit from theta
//       VectorXf thetaPlus(theta);
//       thetaPlus(i) += epsilon;
//       VectorXf thetaMinus(theta);
//       thetaMinus(i) -= epsilon;
//       // Find the objective function value if we add a little to theta
//       VectorXf fvecPlus(values());
//       operator()(thetaPlus, fvecPlus, spde_list,  n_sess, w);
//       // Find the objective function value if we subtract a little from theta
//       VectorXf fvecMinus(values());
//       operator()(thetaMinus, fvecMinus, spde_list,  n_sess, w);
//       // Look at the slope
//       VectorXf fvecDiff(values());
//       fvecDiff = (fvecPlus - fvecMinus) / (2.0f * epsilon);
//
//       fjac.block(0, i, values(), 1) = fvecDiff;
//     }
//
//     return 0;
//   }
//
//   // Number of "data values"
//   int m;
//
//   // Returns 'm', the number of values
//   int values() const { return m; }
//
//   // Number of parameters
//   int n;
//
//   // Returns 'n', the number of inputs
//   int inputs() const { return n; }
// };

// int n = 2;
// int m = 1;
// initFunctor functor_one;
// functor_one.m = m;
// functor_one.n = n;
// VectorXf theta(n);
// LevenbergMarquardt<initFunctor, float> lm(functor_one);
// lm.minimize(theta);

// VectorXf init_theta(VectorXf theta, VectorXf w, List spde_list, int n_sess) {
//   int n = theta.size();
//   int m = 1;
//   initFunctor functor_one;
//   functor_one.m = m;
//   functor_one.n = n;
//   LevenbergMarquardt<initFunctor, float> lm(functor_one);
//   lm.minimize(theta);
//   return theta;
// }
