#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

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
  SimplicialLDLT<Eigen::SparseMatrix<double>> cholQ(Q);
  double lDQ = n_sess * cholQ.vectorD().array().log().sum();
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
};

// [[Rcpp::export]]
Eigen::VectorXd initK(Eigen::VectorXd kappa2, double phi, List in_list, Eigen::VectorXd w, int n_sess) {
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
    VectorXd expKappa2 = kappa2;
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

// [[Rcpp::export]]
double kappa2Obj(double kappa2, List spde, double a_star, double b_star, int n_sess) {
  double lDQ = logDetQt(kappa2, spde, n_sess);
  double out = a_star * kappa2 + b_star / kappa2 - lDQ;
  return out;
}

// [[Rcpp::export]]
double kappa2Brent(double lower, double upper, List spde, double a_star, double b_star, int n_sess, double tol) {
  // Define squared inverse of the golden ratio
  const double c = (3. - sqrt(5.)) / 2.;
  // Initialize local variables
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
  eps = DBL_EPSILON; // machine precision
  tol1 = eps + 1.;
  eps = sqrt(eps);

  a = lower;
  b = upper;
  v = a + c*(b-a);
  x = v;
  x = v;

  d = 0.;
  e = 0.;
  // I don't know what these next three lines mean
  // fx = (*f)(x, info);
  fx = kappa2Obj(x, spde, a_star, b_star, n_sess);
  fv = fx;
  fw = fx;
  tol3 = tol / 3.;

  // Main for loop
  for(;;) {
    xm = (a+b)/2.;
    tol1 = eps * fabs(x) + tol3;
    t2 = tol1 * 2.;
    // Check stopping criterion
    if (fabs(x - xm) <= t2 - (b - a) / 2.) break;
    p = 0.;
    q = 0.;
    r = 0.;
    if (fabs(e) > tol1) { //  fit parabola
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) p = -p; else q = -q;
      r = e;
      e = d;
    }
    if (fabs(p) >= fabs(q * .5 * r) ||
        p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

      if (x < xm) e = b - x; else e = a - x;
      d = c * e;
    }
    else { /* a parabolic-interpolation step */

      d = p / q;
      u = x + d;

      /* f must not be evaluated too close to ax or bx */

      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }

    /* f must not be evaluated too close to x */

    if (fabs(d) >= tol1)
      u = x + d;
    else if (d > 0.)
      u = x + tol1;
    else
      u = x - tol1;

    // fu = (*f)(u, info);
    fu = kappa2Obj(u, spde, a_star, b_star, n_sess);

    /*  update  a, b, v, w, and x */

    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w;    w = x;   x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  // end of main loop
  return x;
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
  int n_spde = w.size();
  n_spde = n_spde / n_sess;
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q = kappa2(0) * Cmat + 2.0 * Gmat + GtCinvG / kappa2(0);
  Eigen::VectorXd Qw = Q * w;
  double wQw = w.transpose() * Qw;
  phi = wQw / (4.0 * M_PI * n_spde * n_sess);
  Eigen::VectorXd theta(2);
  theta(0) = kappa2(0);
  theta(1) = phi;
  Eigen::VectorXd my_diff = theta - old_theta;
  double eps = my_diff.squaredNorm();
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

//' Perform the EM algorithm of the Bayesian GLM fitting
//'
//' @param theta the vector of initial values for theta
//' @param spde a list containing the sparse matrix elements Cmat, Gmat, and GtCinvG
//' @param y the vector of response values
//' @param X the sparse matrix of the data values
//' @param QK a sparse matrix of the prior precision found using the initial values of the hyperparameters
//' @param Psi a sparse matrix representation of the basis function mapping the data locations to the mesh vertices
//' @param A a precomputed matrix crossprod(X%*%Psi)
//' @param Vh A random matrix with elements -1 and 1 used in the Hutchinson estimator of a trace
//' @param tol a value for the tolerance used for a stopping rule (compared to
//'   the squared norm of the differences between theta[s] and theta[s-1])
//' @export
// [[Rcpp::export]]
Rcpp::List findTheta(Eigen::VectorXd theta, List spde, Eigen::VectorXd y, Eigen::SparseMatrix<double> X, Eigen::SparseMatrix<double> QK, Eigen::SparseMatrix<double> Psi, Eigen::SparseMatrix<double> A, Eigen::MatrixXd Vh, double tol) {
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
  // SimplicialLLT<Eigen::SparseMatrix<double>> cholQK;
  // cholQK.compute(QK);
  // cholQK.analyzePattern(QK);
  Eigen::SparseMatrix<double> AdivS2 = A / theta[sig2_ind];
  Eigen::SparseMatrix<double> Sig_inv = QK + AdivS2;
  SimplicialLLT<Eigen::SparseMatrix<double>> cholSigInv;
  cholSigInv.compute(Sig_inv);
  cholSigInv.analyzePattern(Sig_inv);
  Eigen::VectorXd theta_old = theta;
  int Step = 0;
  Rcout << "Step 0 theta: " << theta_old.transpose() << std::endl;
  double eps = tol + 1;
  // Initialize everything
  Eigen::SparseMatrix<double> Xpsi = X * Psi;
  Eigen::VectorXd XpsiY = Xpsi.transpose() * y;
  Eigen::MatrixXd Avh = A * Vh;
  double yy = y.transpose() * y;
  int ySize = y.size();
  double a_star, b_star;
  int idx_start, idx_stop;
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
  double phi_partA, phi_partB, phi_partC;
  double phi_denom = 4.0 * M_PI * n_spde * n_sess;
  double new_kappa2, phi_new;
  while(eps > tol) {
  // while(Step < 1) {
    // cholQK.factorize(QK);
    Eigen::SparseMatrix<double> AdivS2 = A / theta[sig2_ind];
    Eigen::SparseMatrix<double> Sig_inv = QK + AdivS2;
    cholSigInv.factorize(Sig_inv);
    Eigen::VectorXd m = XpsiY / theta(sig2_ind);
    Eigen::VectorXd mu = cholSigInv.solve(m);
    // Solve for sigma_2
    Eigen::VectorXd XpsiMu = Xpsi * mu;
    Eigen::MatrixXd P = cholSigInv.solve(Vh);
    Eigen::MatrixXd PaVh = P.transpose() * Avh;
    Eigen::VectorXd diagPaVh = PaVh.diagonal();
    double TrSigA = diagPaVh.sum() / Ns;
    Eigen::VectorXd Amu = A * mu;
    double muAmu = mu.transpose() * Amu;
    double TrAEww = muAmu + TrSigA;
    double yXpsiMu = y.transpose() * XpsiMu;
    theta[sig2_ind] = (yy - 2 * yXpsiMu + TrAEww) / ySize;
    // Update kappa2 and phi by task
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
        // muCmu
        Cmu = Cmat * muKns;
        muCmu += muKns.transpose() * Cmu;
        // muGmu
        Gmu = Gmat * muKns;
        muGmu += muKns.transpose() * Gmu;
        // muGCGmu
        GCGmu = GtCinvG * muKns;
        muGCGmu += muKns.transpose() * GCGmu;
        // Trace approximations w/ Sigma
        Pkn = P.block(idx_start, 0, n_spde, Ns);
        Vkn = Vh.block(idx_start, 0, n_spde, Ns);
        // Trace of C*Sigma
        CVkn = Cmat * Vkn;
        PCVkn = Pkn.transpose() * CVkn;
        diagPCVkn = PCVkn.diagonal();
        sumDiagPCVkn += diagPCVkn.sum();
        // Trace of G*Sigma
        GVkn = Gmat * Vkn;
        PGVkn = Pkn.transpose() * GVkn;
        diagPGVkn = PGVkn.diagonal();
        sumDiagPGVkn += diagPGVkn.sum();
        // Trace of GCG*Sigma
        GCGVkn = GtCinvG * Vkn;
        PGCGVkn = Pkn.transpose() * GCGVkn;
        diagPGCGVkn = PGCGVkn.diagonal();
        sumDiagPGCGVkn += diagPGCGVkn.sum();
      }
      sumDiagPCVkn = sumDiagPCVkn / Ns;
      sumDiagPGVkn = sumDiagPGVkn / Ns;
      sumDiagPGCGVkn = sumDiagPGCGVkn / Ns;
      // Update kappa2
      a_star = (muCmu + sumDiagPCVkn) / (4.0 * M_PI * theta[k + K]);
      b_star = (muGCGmu + sumDiagPGCGVkn) / (4.0 * M_PI * theta[k + K]);
      // VectorXd new_kappa2 = theta.segment(k,1);
      // new_kappa2 = new_kappa2.log();
      // new_kappa2 = updateK(new_kappa2, spde, a_star, b_star, n_sess);
      // new_kappa2 = new_kappa2.exp();
      new_kappa2 = kappa2Brent(0, 50, spde, a_star, b_star, 2, 1e-4);
      theta[k] = new_kappa2;
      // Update phi
      phi_partA = sumDiagPCVkn + muCmu;
      phi_partA = phi_partA * new_kappa2;
      phi_partB = sumDiagPGVkn + muGmu;
      phi_partB = 2 * phi_partB;
      phi_partC = sumDiagPGCGVkn + muGCGmu;
      phi_partC = phi_partC / new_kappa2;
      double TrQEww = phi_partA + phi_partB + phi_partC;
      phi_new = TrQEww / phi_denom;
      theta[k + K] = phi_new;
      Eigen::SparseMatrix<double> Qk = theta(k) * Cmat + 2.0 * Gmat + GtCinvG / theta(k);
      Qk = Qk / (4.0 * M_PI * theta(k + K));
      for(int ns = 0; ns < n_sess; ns++) {
        int start_i = k * n_spde + ns * K * n_spde;
        setSparseBlock_update(&QK, start_i, start_i, Qk);
      }
    }
    Step += 1;
    Rcout << "Step " << Step << " theta: " << theta.transpose() << std::endl;
    VectorXd thetaDiff = theta - theta_old;
    eps = thetaDiff.squaredNorm();
    theta_old = theta;
  }

  AdivS2 = A / theta[sig2_ind];
  Sig_inv = QK + AdivS2;
  cholSigInv.factorize(Sig_inv);
  Eigen::VectorXd m = XpsiY / theta(sig2_ind);
  Eigen::VectorXd mu = cholSigInv.solve(m);
  List out = List::create(Named("theta_new") = theta,
                          Named("kappa2_new") = theta.segment(0,K),
                          Named("phi_new") = theta.segment(K,K),
                          Named("sigma2_new") = theta(2*K),
                          Named("mu") = mu);
  return out;
}
