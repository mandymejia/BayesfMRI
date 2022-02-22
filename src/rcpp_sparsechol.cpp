#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

//' Find the log of the determinant of Q_tilde
//'
//' @param kappa2 a scalar
//' @param in_list a list with elements Cmat, Gmat, and GtCinvG
//' @param n_sess the integer number of sessions
//' @export
// [[Rcpp::export(rng = false)]]
double logDetQt(double kappa2, const Rcpp::List &in_list, int n_sess) {
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

void makeQt(Eigen::SparseMatrix<double>* Q, double kappa2, const Rcpp::List &spde) {
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  Eigen::SparseMatrix<double> Qt = kappa2 * Cmat + 2. * Gmat + GtCinvG / kappa2;
  for (int k=0; k<Qt.outerSize(); ++k) {
    for (SparseMatrix<double,0,int>::InnerIterator it(Qt,k); it; ++it) {
      Q->coeffRef(it.row(), it.col()) = it.value();
    }
  }
  // return Qt;
}

double kappa2InitObj(double kappa2, double phi, const List &spde, Eigen::VectorXd beta_hat, int n_sess) {
  double lDQ = logDetQt(kappa2, spde, n_sess);
  Eigen::SparseMatrix<double> Cmat = Eigen::SparseMatrix<double> (spde["Cmat"]);
  int n_spde = Cmat.rows();
  Eigen::SparseMatrix<double> Qt(n_spde,n_spde);
  makeQt(&Qt, kappa2, spde);
  Eigen::VectorXd Qw(n_spde), wNs(n_spde);
  double wQw = 0.;
  for(int ns = 0; ns < n_sess; ns++) {
    wNs = beta_hat.segment(ns * n_spde, n_spde);
    Qw = Qt * wNs;
    wQw += wNs.transpose() * Qw;
  }
  wQw = wQw / (4 * M_PI * phi);
  double initObj = wQw - lDQ;
  return initObj;
}

double kappa2BrentInit(double lower, double upper, double phi, const List &spde, Eigen::VectorXd beta_hat, int n_sess, double tol) {
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
  fx = kappa2InitObj(x, phi, spde, beta_hat, n_sess);
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
    fu = kappa2InitObj(u, phi, spde, beta_hat, n_sess);

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

double kappa2Obj(double kappa2, const List &spde, double a_star, double b_star, int n_sess) {
  double lDQ = logDetQt(kappa2, spde, n_sess);
  double out = a_star * kappa2 + b_star / kappa2 - lDQ;
  return out;
}

double kappa2Brent(double lower, double upper, const List &spde, double a_star, double b_star, int n_sess, double tol) {
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
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd initialKP(double kappa2, double phi, List spde, Eigen::VectorXd w, int n_sess, double tol) {
  Eigen::VectorXd old_theta(2);
  old_theta(0) = kappa2;
  old_theta(1) = phi;
  // kappa2 = initK(kappa2, phi, spde, w, n_sess);
  kappa2 = kappa2BrentInit(0., 50., phi, spde, w, n_sess, tol);
  // Rcout << "kappa2: " << kappa2 << std::endl;
  int n_spde = w.size();
  n_spde = n_spde / n_sess;
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q(n_spde,n_spde);
  makeQt(&Q, kappa2, spde);
  // Rcout << "Q dims " << Q.rows() << " X " << Q.cols() << std::endl;
  // Rcout << "w size " << w.size() << std::endl;
  Eigen::VectorXd Qw(n_spde), wNs(n_spde);
  double wQw = 0;
  int start_idx;
  for (int ns = 0; ns < n_sess; ns++) {
    start_idx = ns * n_spde;
    wNs = w.segment(start_idx, n_spde);
    Qw = Q * wNs;
    wQw += wNs.transpose() * Qw;
  }
  phi = wQw / (4.0 * M_PI * n_spde * n_sess);
  Eigen::VectorXd theta(2);
  theta(0) = kappa2;
  theta(1) = phi;
  Eigen::VectorXd my_diff = theta - old_theta;
  double eps = my_diff.squaredNorm();
  int num_steps = 1;
  while(eps > tol) {
    old_theta(0) = kappa2;
    old_theta(1) = phi;
    wQw = 0.;
    kappa2 = kappa2BrentInit(0., 50., phi, spde, w, n_sess, tol);
    makeQt(&Q, kappa2, spde);
    for (int ns = 0; ns < n_sess; ns++) {
      start_idx = ns * n_spde;
      wNs = w.segment(start_idx, n_spde);
      Qw = Q * wNs;
      wQw += wNs.transpose() * Qw;
    }
    phi = wQw / (4.0 * M_PI * n_spde * n_sess);
    theta(0) = kappa2;
    theta(1) = phi;
    my_diff = theta - old_theta;
    eps = my_diff.squaredNorm();
    num_steps += 1;
  }
  return theta;
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

Eigen::VectorXd theta_fixpt(Eigen::VectorXd theta, const Eigen::SparseMatrix<double> A,
                            Eigen::SparseMatrix<double> QK, SimplicialLLT<Eigen::SparseMatrix<double>> &cholSigInv,
                            const Eigen::VectorXd XpsiY, const Eigen::SparseMatrix<double> Xpsi,
                            const Eigen::MatrixXd Vh, const Eigen::MatrixXd Avh,
                            const Eigen::VectorXd y, const double yy, const List spde) {
  // Bring in the spde matrices
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  // Grab metadata
  int K = (theta.size() - 1) / 2;
  int sig2_ind = theta.size() - 1;
  int nKs = A.rows();
  int ySize = y.size();
  int n_spde = Cmat.rows();
  int n_sess = nKs / (n_spde * K);
  int Ns = Vh.cols();
  // Initialize objects
  Eigen::SparseMatrix<double> AdivS2(nKs,nKs), Sig_inv(nKs,nKs), Qk(n_spde, n_spde);
  Eigen::VectorXd theta_new = theta;
  Eigen::VectorXd muKns(n_spde), Cmu(n_spde), Gmu(n_spde), diagPCVkn(Ns);
  Eigen::VectorXd diagPGVkn(Ns), GCGmu(n_spde), diagPGCGVkn(Ns);
  Eigen::MatrixXd Pkn(n_spde, Ns), Vkn(n_spde, Ns), CVkn(n_spde, Ns), PCVkn(Ns, Ns);
  Eigen::MatrixXd GVkn(n_spde, Ns), PGVkn(Ns, Ns), GCGVkn(n_spde, Ns), PGCGVkn(Ns,Ns);
  double a_star, b_star, muCmu, muGmu, sumDiagPCVkn, sumDiagPGVkn, muGCGmu;
  double sumDiagPGCGVkn, phi_partA, phi_partB, phi_partC, new_kappa2, phi_new;
  double phi_denom = 4.0 * M_PI * n_spde * n_sess;
  int idx_start;
  // Begin update
  for(int k = 0; k < K ; k++) {
    makeQt(&Qk, theta(k), spde);
    // Qk = theta(k) * Cmat + 2.0 * Gmat + GtCinvG / theta(k);
    Qk = Qk / (4.0 * M_PI * theta(k + K));
    for(int ns = 0; ns < n_sess; ns++) {
      int start_i = k * n_spde + ns * K * n_spde;
      setSparseBlock_update(&QK, start_i, start_i, Qk);
    }
  }
  AdivS2 = A / theta[sig2_ind];
  Sig_inv = QK + AdivS2;
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
  theta_new[sig2_ind] = (yy - 2 * yXpsiMu + TrAEww) / ySize;
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
      // idx_stop = idx_start + n_spde;
      muKns = mu.segment(idx_start,n_spde);
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
    new_kappa2 = kappa2Brent(0, 50, spde, a_star, b_star, 2, 1e-4);
    theta_new[k] = new_kappa2;
    // Update phi
    phi_partA = sumDiagPCVkn + muCmu;
    phi_partA = phi_partA * new_kappa2;
    phi_partB = sumDiagPGVkn + muGmu;
    phi_partB = 2 * phi_partB;
    phi_partC = sumDiagPGCGVkn + muGCGmu;
    phi_partC = phi_partC / new_kappa2;
    double TrQEww = phi_partA + phi_partB + phi_partC;
    phi_new = TrQEww / phi_denom;
    theta_new[k + K] = phi_new;
  }
  return(theta_new);
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
// [[Rcpp::export(rng = false)]]
Rcpp::List findTheta(Eigen::VectorXd theta, List spde, Eigen::VectorXd y,
                     Eigen::SparseMatrix<double> X, Eigen::SparseMatrix<double> QK,
                     Eigen::SparseMatrix<double> Psi, Eigen::SparseMatrix<double> A,
                     Eigen::MatrixXd Vh, double tol) {
  // Bring in the spde matrices
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  // int Ns = Vh.cols(); // For the Hutchinson estimator involving the trace of matrix products with Sigma
  int n_spde = Cmat.rows();
  // int nKs = A.rows();
  int K = theta.size();
  K = (K - 1) / 2;
  // int n_sess = nKs / (n_spde * K);
  int sig2_ind = 2*K;
  Eigen::SparseMatrix<double> AdivS2 = A / theta[sig2_ind];
  Eigen::SparseMatrix<double> Sig_inv = QK + AdivS2;
  SimplicialLLT<Eigen::SparseMatrix<double>> cholSigInv;
  cholSigInv.compute(Sig_inv);
  cholSigInv.analyzePattern(Sig_inv);
  // Eigen::VectorXd theta_old = theta;
  int Step = 0;
  Rcout << "Step 0 theta: " << theta.transpose() << std::endl;
  double eps = tol + 1;
  // Initialize everything
  Eigen::SparseMatrix<double> Xpsi = X * Psi;
  Eigen::SparseMatrix<double> Qk(n_spde, n_spde);
  Eigen::VectorXd XpsiY = Xpsi.transpose() * y;
  Eigen::MatrixXd Avh = A * Vh;
  double yy = y.transpose() * y;
  Eigen::VectorXd theta_new;
  // int ySize = y.size();
  // double a_star, b_star, muCmu, muGmu, sumDiagPCVkn, sumDiagPGVkn, muGCGmu;
  // double sumDiagPGCGVkn, phi_partA, phi_partB, phi_partC, new_kappa2, phi_new;
  // double phi_denom = 4.0 * M_PI * n_spde * n_sess;
  // double TrSigA, muAmu, TrAEww, yXpsiMu;
  // int idx_start;
  // Eigen::VectorXd muKns(n_spde), Cmu(n_spde), Gmu(n_spde), diagPCVkn(Ns);
  // Eigen::VectorXd diagPGVkn(Ns), GCGmu(n_spde), diagPGCGVkn(Ns), thetaDiff((2*K + 1));
  // Eigen::VectorXd m(nKs), mu(nKs), XpsiMu(ySize), diagPaVh(Ns), Amu(nKs);
  // Eigen::MatrixXd Pkn(n_spde, Ns), Vkn(n_spde, Ns), CVkn(n_spde, Ns), PCVkn(Ns, Ns);
  // Eigen::MatrixXd GVkn(n_spde, Ns), PGVkn(Ns, Ns), GCGVkn(n_spde, Ns), PGCGVkn(Ns,Ns);
  // Eigen::MatrixXd P(nKs, Ns), PaVh(Ns,Ns);
  Eigen::VectorXd thetaDiff(2 * K + 1);
  while(eps > tol) {
    theta_new = theta_fixpt(theta, A, QK, cholSigInv, XpsiY, Xpsi, Vh, Avh, y, yy, spde);
  // // while(Step < 1) {
  //   AdivS2 = A / theta[sig2_ind];
  //   Sig_inv = QK + AdivS2;
  //   cholSigInv.factorize(Sig_inv);
  //   m = XpsiY / theta(sig2_ind);
  //   mu = cholSigInv.solve(m);
  //   // Solve for sigma_2
  //   XpsiMu = Xpsi * mu;
  //   P = cholSigInv.solve(Vh);
  //   PaVh = P.transpose() * Avh;
  //   diagPaVh = PaVh.diagonal();
  //   TrSigA = diagPaVh.sum() / Ns;
  //   Amu = A * mu;
  //   muAmu = mu.transpose() * Amu;
  //   TrAEww = muAmu + TrSigA;
  //   yXpsiMu = y.transpose() * XpsiMu;
  //   theta[sig2_ind] = (yy - 2 * yXpsiMu + TrAEww) / ySize;
  //   // Update kappa2 and phi by task
  //   for(int k = 0; k < K; k++) {
  //     a_star = 0.0;
  //     b_star = 0.0;
  //     muCmu = 0.0;
  //     muGmu = 0.0;
  //     muGCGmu = 0.0;
  //     sumDiagPCVkn = 0.0;
  //     sumDiagPGVkn = 0.0;
  //     sumDiagPGCGVkn = 0.0;
  //     for(int ns = 0; ns < n_sess; ns++) {
  //       idx_start = k * n_spde + ns * K * n_spde;
  //       // idx_stop = idx_start + n_spde;
  //       muKns = mu.segment(idx_start,n_spde);
  //       // muCmu
  //       Cmu = Cmat * muKns;
  //       muCmu += muKns.transpose() * Cmu;
  //       // muGmu
  //       Gmu = Gmat * muKns;
  //       muGmu += muKns.transpose() * Gmu;
  //       // muGCGmu
  //       GCGmu = GtCinvG * muKns;
  //       muGCGmu += muKns.transpose() * GCGmu;
  //       // Trace approximations w/ Sigma
  //       Pkn = P.block(idx_start, 0, n_spde, Ns);
  //       Vkn = Vh.block(idx_start, 0, n_spde, Ns);
  //       // Trace of C*Sigma
  //       CVkn = Cmat * Vkn;
  //       PCVkn = Pkn.transpose() * CVkn;
  //       diagPCVkn = PCVkn.diagonal();
  //       sumDiagPCVkn += diagPCVkn.sum();
  //       // Trace of G*Sigma
  //       GVkn = Gmat * Vkn;
  //       PGVkn = Pkn.transpose() * GVkn;
  //       diagPGVkn = PGVkn.diagonal();
  //       sumDiagPGVkn += diagPGVkn.sum();
  //       // Trace of GCG*Sigma
  //       GCGVkn = GtCinvG * Vkn;
  //       PGCGVkn = Pkn.transpose() * GCGVkn;
  //       diagPGCGVkn = PGCGVkn.diagonal();
  //       sumDiagPGCGVkn += diagPGCGVkn.sum();
  //     }
  //     sumDiagPCVkn = sumDiagPCVkn / Ns;
  //     sumDiagPGVkn = sumDiagPGVkn / Ns;
  //     sumDiagPGCGVkn = sumDiagPGCGVkn / Ns;
  //     // Update kappa2
  //     a_star = (muCmu + sumDiagPCVkn) / (4.0 * M_PI * theta[k + K]);
  //     b_star = (muGCGmu + sumDiagPGCGVkn) / (4.0 * M_PI * theta[k + K]);
  //     new_kappa2 = kappa2Brent(0, 50, spde, a_star, b_star, 2, 1e-4);
  //     theta[k] = new_kappa2;
  //     // Update phi
  //     phi_partA = sumDiagPCVkn + muCmu;
  //     phi_partA = phi_partA * new_kappa2;
  //     phi_partB = sumDiagPGVkn + muGmu;
  //     phi_partB = 2 * phi_partB;
  //     phi_partC = sumDiagPGCGVkn + muGCGmu;
  //     phi_partC = phi_partC / new_kappa2;
  //     double TrQEww = phi_partA + phi_partB + phi_partC;
  //     phi_new = TrQEww / phi_denom;
  //     theta[k + K] = phi_new;
  //     makeQt(&Qk, theta(k), spde);
  //     // Qk = theta(k) * Cmat + 2.0 * Gmat + GtCinvG / theta(k);
  //     Qk = Qk / (4.0 * M_PI * theta(k + K));
  //     for(int ns = 0; ns < n_sess; ns++) {
  //       int start_i = k * n_spde + ns * K * n_spde;
  //       setSparseBlock_update(&QK, start_i, start_i, Qk);
  //     }
  //   }
    Step += 1;
    Rcout << "Step " << Step << " theta: " << theta_new.transpose() << std::endl;
    thetaDiff = theta_new - theta;
    eps = thetaDiff.squaredNorm();
    theta = theta_new;
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

