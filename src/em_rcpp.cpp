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
//' 
// [[Rcpp::export(.logDetQt, rng = false)]]
double logDetQt(double kappa2, const Rcpp::List &in_list, double n_sess) {
  // Load parameters
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (in_list["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (in_list["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (in_list["GtCinvG"]);
  // Create SparseMatrix Q
  Eigen::SparseMatrix<double> Q= kappa2 * Cmat + 2.0 * Gmat + GtCinvG / kappa2;
  SimplicialLDLT<Eigen::SparseMatrix<double> > cholQ(Q);
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

double kappa2InitObj(double kappa2, double phi, const List &spde, Eigen::VectorXd beta_hat, double n_sess) {
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
  // Rcout << "log_det_Q = " << lDQ << ", bQb = " << wQw << std::endl;
  wQw = wQw / (8. * M_PI * phi);
  lDQ = lDQ / 2.;
  double initObj = wQw - lDQ;
  return initObj;
}

double kappa2BrentInit(double lower, double upper, double phi, const List &spde, Eigen::VectorXd beta_hat, double n_sess) {
  // Define squared inverse of the golden ratio
  const double c = (3. - std::sqrt(5.)) / 2.;
  // Initialize local variables
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
  eps = DBL_EPSILON; // machine precision
  tol1 = eps + 1.;
  eps = std::sqrt(eps);
  double tol = std::sqrt(eps);

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
    tol1 = eps * std::abs(x) + tol3;
    t2 = tol1 * 2.;
    // Check stopping criterion
    if (std::abs(x - xm) <= t2 - (b - a) / 2.) break;
    p = 0.;
    q = 0.;
    r = 0.;
    if (std::abs(e) > tol1) { //  fit parabola
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) p = -p; else q = -q;
      r = e;
      e = d;
    }
    if (std::abs(p) >= std::abs(q * .5 * r) ||
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

    if (std::abs(d) >= tol1)
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
  // Rcout << "objective = " << fx << std::endl;
  return x;
}

double kappa2Obj(double kappa2, const List &spde, double a_star, double b_star, double n_sess) {
  double lDQ = logDetQt(kappa2, spde, n_sess);
  double out = a_star * kappa2 + b_star / kappa2 - lDQ;
  return out;
}

double kappa2Brent(double lower, double upper, const List &spde, double a_star, double b_star, double n_sess) {
  // Define squared inverse of the golden ratio
  const double c = (3. - std::sqrt(5.)) / 2.;
  // Initialize local variables
  double a, b, d, e, p, q, r, u, v, w, x;
  double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
  eps = DBL_EPSILON; // machine precision
  tol1 = eps + 1.;
  eps = std::sqrt(eps);
  double tol = std::sqrt(eps);

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
    tol1 = eps * std::abs(x) + tol3;
    t2 = tol1 * 2.;
    // Check stopping criterion
    if (std::abs(x - xm) <= t2 - (b - a) / 2.) break;
    p = 0.;
    q = 0.;
    r = 0.;
    if (std::abs(e) > tol1) { //  fit parabola
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.;
      if (q > 0.) p = -p; else q = -q;
      r = e;
      e = d;
    }
    if (std::abs(p) >= std::abs(q * .5 * r) ||
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

    if (std::abs(d) >= tol1)
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
  // Rcout << "objective = " << fx;
  return x;
}

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

Eigen::VectorXd init_fixptC(Eigen::VectorXd theta, Eigen::VectorXd w, List spde, double n_sess) {
  int n_spde = w.size();
  int start_idx;
  Eigen::VectorXd wNs(n_spde);
  n_spde = n_spde / n_sess;
  Eigen::VectorXd Qw(n_spde);
  double wQw = 0.;
  // theta(0) = kappa2BrentInit(0., 50., theta(1), spde, w, n_sess, tol);
  theta(0) = kappa2BrentInit(0., 50., theta(1), spde, w, n_sess);
  Eigen::SparseMatrix<double> Q(n_spde,n_spde);
  makeQt(&Q, theta(0), spde);
  for (int ns = 0; ns < n_sess; ns++) {
    start_idx = ns * n_spde;
    wNs = w.segment(start_idx, n_spde);
    Qw = Q * wNs;
    wQw += wNs.transpose() * Qw;
  }
  theta(1) = wQw / (4.0 * M_PI * n_spde * n_sess);
  return theta;
}

SquaremOutput init_squarem2(Eigen::VectorXd par, Eigen::VectorXd w, List spde, double n_sess, double tol){
  double res,parnorm,kres;
  Eigen::VectorXd pcpp,p1cpp,p2cpp,pnew,ptmp;
  Eigen::VectorXd q1,q2,sr2,sq2,sv2,srv;
  Eigen::VectorXd diffp1p, diffp2p1;
  double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;
  int iter,feval;
  bool conv,extrap;
  stepmin=SquaremDefault.stepmin0;
  stepmax=SquaremDefault.stepmax0;
  if(SquaremDefault.trace){Rcout<<"Squarem-2"<<std::endl;}

  iter=1;pcpp=par;pnew=par;
  feval=0;conv=true;

  const long int parvectorlength=pcpp.size();

  while(feval<SquaremDefault.maxiter){
    //Step 1
    extrap = true;
    // try{p1cpp=fixptfn(pcpp);feval++;}
    try{p1cpp=init_fixptC(pcpp, w, spde, n_sess);feval++;}
    catch(...){
      Rcout<<"Error in fixptfn function evaluation";
      return sqobjnull;
    }

    // diffp1p = p1cpp - pcpp;
    // sr2_scalar = diffpp1p.squaredNorm();
    sr2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sr2_scalar+=std::pow(p1cpp[i]-pcpp[i],2.0);}
    if(std::sqrt(sr2_scalar)<SquaremDefault.tol){break;}

    //Step 2
    try{p2cpp=init_fixptC(p1cpp,  w, spde, n_sess);feval++;}
    catch(...){
      Rcout<<"Error in fixptfn function evaluation";
      return sqobjnull;
    }
    // diffp2p1 = p2cpp - p1cpp;
    // sq2_scalar= diffp2p1.squaredNorm();
    sq2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sq2_scalar+=std::pow(p2cpp[i]-p1cpp[i],2.0);}
    sq2_scalar=std::sqrt(sq2_scalar);
    if (sq2_scalar<SquaremDefault.tol){break;}
    res=sq2_scalar;

    // p2M2p1Pp = p2cpp - 2. * p1cpp + pcpp;
    // sv2_scalar = p2M2p1Pp.squaredNorm();
    sv2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sv2_scalar+=std::pow(p2cpp[i]-2*p1cpp[i]+pcpp[i],2.0);}
    // p1Mp = p1cpp - pcpp;
    // p2M2p1PpTp1Mp = p2M2p1Pp * p1Mp;
    // srv_scalar = p2M2p1PpTp1Mp.squaredNorm();
    srv_scalar=0;
    for (int i=0;i<parvectorlength;i++){srv_scalar+=(p2cpp[i]-2*p1cpp[i]+pcpp[i])*(p1cpp[i]-pcpp[i]);}
    //std::cout<<"sr2,sv2,srv="<<sr2_scalar<<","<<sv2_scalar<<","<<srv_scalar<<std::endl;//debugging

    //Step 3 Proposing new value
    switch (SquaremDefault.method){
    case 1: alpha= -srv_scalar/sv2_scalar;
    case 2: alpha= -sr2_scalar/srv_scalar;
    case 3: alpha= std::sqrt(sr2_scalar/sv2_scalar);
    }

    alpha=std::max(stepmin,std::min(stepmax,alpha));
    //std::cout<<"alpha="<<alpha<<std::endl;//debugging
    for (int i=0;i<parvectorlength;i++){pnew[i]=pcpp[i]+2.0*alpha*(p1cpp[i]-pcpp[i])+alpha*alpha*(p2cpp[i]-2.0*p1cpp[i]+pcpp[i]);}
    // pnew = pcpp + 2.0*alpha*q1 + alpha*alpha*(q2-q1);

    //Step 4 stabilization
    if(std::abs(alpha-1)>0.01){
      try{ptmp=init_fixptC(pnew,  w, spde, n_sess);feval++;}
      catch(...){
        pnew=p2cpp;
        if(alpha==stepmax){
          stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
        }
        alpha=1;
        extrap=false;
        if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
        if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}
        pcpp=pnew;
        if(SquaremDefault.trace){Rcout<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
        iter++;
        continue;//next round in while loop
      }
      res=0;
      for (int i=0;i<parvectorlength;i++){res+=std::pow(ptmp[i]-pnew[i],2.0);}
      res=std::sqrt(res);
      parnorm=0;
      for (int i=0;i<parvectorlength;i++){parnorm+=std::pow(p2cpp[i],2.0);}
      parnorm=std::sqrt(parnorm/parvectorlength);
      kres=SquaremDefault.kr*(1+parnorm)+sq2_scalar;
      if(res <= kres){
        pnew=ptmp;
      }else{
        pnew=p2cpp;
        if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
        alpha=1;
        extrap=false;
      }
    }

    if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
    if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}

    pcpp=pnew;
    if(SquaremDefault.trace){Rcout<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
    iter++;
  }

  if (feval >= SquaremDefault.maxiter){conv=false;}

  //assigning values
  sqobj.par=pcpp;
  sqobj.valueobjfn=NAN;
  sqobj.iter=iter;
  sqobj.pfevals=feval;
  sqobj.objfevals=0;
  sqobj.convergence=conv;
  return(sqobj);
}

//' Find the initial values of kappa2 and phi
//'
//' @param theta a vector of length two containing the range and scale parameters
//'   kappa2 and phi, in that order
//' @param spde a list containing the sparse matrix elements Cmat, Gmat, and GtCinvG
//' @param w the beta_hat estimates for a single task
//' @param n_sess the number of sessions
//' @param tol the stopping rule tolerance
//' @param verbose (logical) Should intermediate output be displayed?
//' 
// [[Rcpp::export(.initialKP, rng = false)]]
Eigen::VectorXd initialKP(Eigen::VectorXd theta, List spde, Eigen::VectorXd w,
                          double n_sess, double tol, bool verbose) {
  int n_spde = w.size();
  n_spde = n_spde / n_sess;
  // Set up implementation without EM
  // double eps = tol + 1;
  // Eigen::VectorXd new_theta(theta.size()), diffTheta(theta.size());
  // while(eps > tol) {
  //   new_theta = init_fixptC(theta, w, spde, n_sess);
  //   diffTheta = new_theta - theta;
  //   eps = diffTheta.squaredNorm();
  //   eps = sqrt(eps);
  //   theta = new_theta;
  // }
  // Implementation with in EM
  SquaremOutput SQ_out;
  SquaremDefault.tol = tol;
  SquaremDefault.trace = verbose;
  SQ_out = init_squarem2(theta, w, spde, n_sess, tol);
  // Rcout << "valueobjfn = " << SQ_out.valueobjfn << ", iter = " << SQ_out.iter;
  // Rcout << ", fpevals = " << SQ_out.pfevals << ", objevals = " << SQ_out.objfevals;
  // Rcout << ", convergence = " << SQ_out.convergence << std::endl;
  theta= SQ_out.par;
  return theta;
}

Eigen::MatrixXd makeV(int n_spde, int Ns) {
  Eigen::MatrixXd V(n_spde,Ns);
  Rcpp::NumericVector x(1);
  for(int i = 0; i < n_spde; i++) {
    for(int j = 0; j < Ns; j++) {
      // x = std::rand() / ((RAND_MAX + 1u) / 2);
      x = Rcpp::runif(1);
      if(x[0] <= 0.5) {
        V(i,j) = -1;
      } else {
        V(i,j) = 1;
      }
    }
  }
  return V;
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

double emObj(Eigen::VectorXd theta, const Eigen::SparseMatrix<double> A,
             Eigen::SparseMatrix<double> QK,
             SimplicialLLT<Eigen::SparseMatrix<double> > &cholSigInv,
             const Eigen::VectorXd XpsiY, const Eigen::SparseMatrix<double> Xpsi,
             const int Ns, const Eigen::VectorXd y, const List spde) {
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
  double n_sess = nKs / (n_spde * K);
  Eigen::MatrixXd Vh = makeV(nKs,Ns);
  // Initialize objects
  Eigen::SparseMatrix<double> AdivS2(nKs,nKs), Sig_inv(nKs,nKs), Qk(n_spde, n_spde);
  for(int k = 0; k < K ; k++) {
    makeQt(&Qk, theta(k), spde);
    Qk = Qk / (4.0 * M_PI * theta(k + K));
    for(int ns = 0; ns < n_sess; ns++) {
      int start_i = k * n_spde + ns * K * n_spde;
      setSparseBlock_update(&QK, start_i, start_i, Qk);
    }
  }
  AdivS2 = A / theta[sig2_ind];
  Sig_inv = QK + AdivS2;
  cholSigInv.factorize(Sig_inv);
  Eigen::MatrixXd P = cholSigInv.solve(Vh);
  Eigen::MatrixXd PqVh = P.transpose() * QK;
  Eigen::VectorXd diagPqVh = PqVh.diagonal();
  double TrSigQ = diagPqVh.sum() / Ns;
  Eigen::VectorXd m = XpsiY / theta(sig2_ind);
  Eigen::VectorXd mu = cholSigInv.solve(m);
  Eigen::MatrixXd muTmu = mu.transpose() * mu;
  Eigen::MatrixXd QmuTmu = QK * muTmu;
  double TrQmuTmu = QmuTmu.diagonal().sum();
  SimplicialLDLT<Eigen::SparseMatrix<double> > cholQ(QK);
  double lDQ = cholQ.vectorD().array().log().sum();
  Eigen::VectorXd XB = Xpsi * mu;
  Eigen::VectorXd y_til = y - XB;
  double ytil2 = y_til.transpose() * y_til;
  double llik = -ySize * log(theta[sig2_ind]);
  llik = llik - ytil2 / theta[sig2_ind];
  llik = llik - lDQ;
  llik = llik - TrSigQ - TrQmuTmu;
  llik = llik / 2.;
  return llik;
}

Eigen::VectorXd theta_fixpt(Eigen::VectorXd theta, const Eigen::SparseMatrix<double> A,
                            Eigen::SparseMatrix<double> QK, SimplicialLLT<Eigen::SparseMatrix<double> > &cholSigInv,
                            const Eigen::VectorXd XpsiY, const Eigen::SparseMatrix<double> Xpsi,
                            const int Ns, const Eigen::VectorXd y,
                            const double yy, const List spde, double tol) {
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
  double n_sess = nKs / (n_spde * K);
  // int Ns = Vh.cols();
  Eigen::MatrixXd Vh = makeV(nKs,Ns);
  // Rcout << "dim(A)" << A.rows() << " x " << A.cols() << ", dim(Vh) = " << Vh.rows() << " x " << Vh.cols() << std::endl;
  Eigen::MatrixXd Avh = A * Vh;
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
  // Rcout << "First 6 values of mu: " << mu.segment(0,6).transpose() << std::endl;
  // Solve for sigma_2
  Eigen::VectorXd XpsiMu = Xpsi * mu;
  Eigen::MatrixXd P = cholSigInv.solve(Vh);
  Eigen::MatrixXd PaVh = P.transpose() * Avh;
  Eigen::VectorXd diagPaVh = PaVh.diagonal();
  double TrSigA = diagPaVh.sum() / Ns;
  Eigen::VectorXd Amu = A * mu;
  double muAmu = mu.transpose() * Amu;
  double TrAEww = muAmu + TrSigA;
  // Rcout << "TrAEww = " << TrAEww << std::endl;
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
    // Rcout << "muCmu = " << muCmu << ", muGCGmu = " << muGCGmu;
    // Rcout << ", TrCSig = " << sumDiagPCVkn << ", TrGCGSig = " << sumDiagPGCGVkn << std::endl;
    // Rcout << "head(diagPCV) = " << diagPCVkn.segment(0,6).transpose() << ", head(diagPGCGV) = " << diagPGCGVkn.segment(0,6).transpose() << std::endl;
    a_star = (muCmu + sumDiagPCVkn) / (4.0 * M_PI * theta[k + K]);
    b_star = (muGCGmu + sumDiagPGCGVkn) / (4.0 * M_PI * theta[k + K]);
    // Rcout << "k = " << k << " a_star = " << a_star << " b_star = " << b_star << std::endl;
    new_kappa2 = kappa2Brent(0., 50., spde, a_star, b_star, n_sess);
    // Rcout << ", new_kappa2 = " << new_kappa2 << std::endl;
    theta_new[k] = new_kappa2;
    // Update phi
    phi_partA = sumDiagPCVkn + muCmu;
    phi_partA = phi_partA * new_kappa2;
    phi_partB = sumDiagPGVkn + muGmu;
    phi_partB = 2 * phi_partB;
    phi_partC = sumDiagPGCGVkn + muGCGmu;
    phi_partC = phi_partC / new_kappa2;
    double TrQEww = phi_partA + phi_partB + phi_partC;
    // Rcout << "TrQEww = " << TrQEww << std::endl;
    phi_new = TrQEww / phi_denom;
    theta_new[k + K] = phi_new;
  }
  return(theta_new);
}

// #include <iostream>
// #include <algorithm>
// #include <cmath>
// #include <math.h>
// #include <vector>
// #include <numeric>

// using namespace std;



SquaremOutput theta_squarem2(Eigen::VectorXd par, const Eigen::SparseMatrix<double> A,
                       Eigen::SparseMatrix<double> QK, SimplicialLLT<Eigen::SparseMatrix<double> > &cholSigInv,
                       const Eigen::VectorXd XpsiY, const Eigen::SparseMatrix<double> Xpsi,
                       const int Ns, const Eigen::VectorXd y, const double yy,
                       const List spde, double tol, bool verbose){
  double res,parnorm,kres;;//, theta_length=par.size(); //unused
  Eigen::VectorXd pcpp,p1cpp,p2cpp,pnew,ptmp;
  Eigen::VectorXd q1,q2,sr2,sq2,sv2,srv;
  double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;
  // double ob_pcpp, ob_p1cpp, ob_p2cpp, ob_ptmp, ob_pnew, rel_llik_pp1, rel_llik_p1p2, rel_llik_tmpNew;
  int iter,feval;
  bool conv,extrap;
  stepmin=SquaremDefault.stepmin0;
  stepmax=SquaremDefault.stepmax0;
  if(verbose){Rcout<<"Squarem-2"<<std::endl;}

  iter=1;pcpp=par;pnew=par;
  feval=0;conv=true;
  // ob_pcpp = emObj(pcpp,A,QK,cholSigInv,XpsiY,Xpsi,Ns,y,spde);

  const long int parvectorlength=pcpp.size();

  while(feval<SquaremDefault.maxiter){
    //Step 1
    extrap = true;
    // try{p1cpp=fixptfn(pcpp);feval++;}
    try{p1cpp=theta_fixpt(pcpp, A, QK, cholSigInv, XpsiY, Xpsi, Ns, y, yy, spde, tol);feval++;}
    catch(...){
      Rcout<<"Error in fixptfn function evaluation";
      return sqobjnull;
    }
    // ob_p1cpp = emObj(p1cpp,A,QK,cholSigInv,XpsiY,Xpsi,Ns,y,spde);
    // rel_llik_pp1 = std::abs(ob_p1cpp - ob_pcpp) / std::abs(ob_pcpp);

    sr2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sr2_scalar+=std::pow(p1cpp[i]-pcpp[i],2.0);}
    // sr2_scalar = std::sqrt(sr2_scalar);
    double pcpp_norm = pcpp.squaredNorm();
    pcpp_norm = std::sqrt(pcpp_norm);
    // if(sr2_scalar / pcpp_norm < tol){break;}
    // if(std::sqrt(sr2_scalar)<tol){break;}
    if(std::sqrt(sr2_scalar) / pcpp_norm <tol){break;}
    // if(rel_llik_pp1<tol){break;}

    //Step 2
    try{p2cpp=theta_fixpt(p1cpp, A, QK, cholSigInv, XpsiY, Xpsi, Ns, y, yy, spde, tol);feval++;}
    catch(...){
      Rcout<<"Error in fixptfn function evaluation";
      return sqobjnull;
    }
    // ob_p2cpp = emObj(p2cpp,A,QK,cholSigInv,XpsiY,Xpsi,Ns,y,spde);
    // rel_llik_p1p2 = std::abs(ob_p2cpp - ob_p1cpp) / std::abs(ob_p1cpp);
    sq2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sq2_scalar+=std::pow(p2cpp[i]-p1cpp[i],2.0);}
    sq2_scalar=std::sqrt(sq2_scalar);
    double p1cpp_norm = p1cpp.squaredNorm();
    p1cpp_norm = std::sqrt(p1cpp_norm);
    if (sq2_scalar / p1cpp_norm <tol){break;}
    // if (rel_llik_p1p2<tol){break;}
    // if (sq2_scalar<tol){break;}
    res=sq2_scalar;
    // res=rel_llik_p1p2;

    sv2_scalar=0;
    for (int i=0;i<parvectorlength;i++){sv2_scalar+=std::pow(p2cpp[i]-2.0*p1cpp[i]+pcpp[i],2.0);}
    srv_scalar=0;
    for (int i=0;i<parvectorlength;i++){srv_scalar+=(p2cpp[i]-2.0*p1cpp[i]+pcpp[i])*(p1cpp[i]-pcpp[i]);}
    //std::cout<<"sr2,sv2,srv="<<sr2_scalar<<","<<sv2_scalar<<","<<srv_scalar<<std::endl;//debugging

    //Step 3 Proposing new value
    switch (SquaremDefault.method){
    case 1: alpha= -srv_scalar/sv2_scalar;
    case 2: alpha= -sr2_scalar/srv_scalar;
    case 3: alpha= std::sqrt(sr2_scalar/sv2_scalar);
    }

    alpha=std::max(stepmin,std::min(stepmax,alpha));
    //std::cout<<"alpha="<<alpha<<std::endl;//debugging
    for (int i=0;i<parvectorlength;i++){pnew[i]=pcpp[i]+2.0*alpha*(p1cpp[i]-pcpp[i])+alpha*alpha*(p2cpp[i]-2.0*p1cpp[i]+pcpp[i]);}
    // pnew = pcpp + 2.0*alpha*q1 + alpha*alpha*(q2-q1);
    // ob_pnew = emObj(pnew,A,QK,cholSigInv,XpsiY,Xpsi,Ns,y,spde);

    //Step 4 stabilization
    if(std::abs(alpha-1)>0.01){
      try{ptmp=theta_fixpt(pnew, A, QK, cholSigInv, XpsiY, Xpsi, Ns, y, yy, spde, tol);feval++;}
      catch(...){
        pnew=p2cpp;
        if(alpha==stepmax){
          stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
        }
        alpha=1;
        extrap=false;
        if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
        if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}
        pcpp=pnew;
        if(verbose){Rcout<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
        iter++;
        continue;//next round in while loop
      }
      // ob_ptmp = emObj(ptmp,A,QK,cholSigInv,XpsiY,Xpsi,Ns,y,spde);
      res=0;
      for (int i=0;i<parvectorlength;i++){res+=std::pow(ptmp[i]-pnew[i],2.0);}
      res=std::sqrt(res);
      // double pnew_norm = pnew.squaredNorm();
      // pnew_norm = std::sqrt(pnew_norm);
      // res = res/ pnew_norm;
      // rel_llik_tmpNew = std::abs(ob_ptmp - ob_pnew)/ std::abs(ob_pnew);
      parnorm=0;
      for (int i=0;i<parvectorlength;i++){parnorm+=std::pow(p2cpp[i],2.0);}
      parnorm=std::sqrt(parnorm/parvectorlength);
      kres=SquaremDefault.kr*(1+parnorm)+sq2_scalar;
      if(res <= kres){
        pnew=ptmp;
      }else{
        pnew=p2cpp;
        if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
        alpha=1;
        extrap=false;
      }
      // res = rel_llik_tmpNew;
    }

    if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
    if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}

    pcpp=pnew;
    if(verbose){Rcout<<"Residual: "<<res<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
    iter++;
  }

  if (feval >= SquaremDefault.maxiter){conv=false;}

  //assigning values
  sqobj.par=pcpp;
  sqobj.valueobjfn=NAN;
  sqobj.iter=iter;
  sqobj.pfevals=feval;
  sqobj.objfevals=0;
  sqobj.convergence=conv;
  return(sqobj);
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
//' @param Ns the number of columns for the random matrix used in the Hutchinson estimator
//' @param tol a value for the tolerance used for a stopping rule (compared to
//'   the squared norm of the differences between \code{theta(s)} and \code{theta(s-1)})
//' @param verbose (logical) Should intermediate output be displayed?
//' 
// [[Rcpp::export(.findTheta, rng = false)]]
Rcpp::List findTheta(Eigen::VectorXd theta, List spde, Eigen::VectorXd y,
                     Eigen::SparseMatrix<double> X, Eigen::SparseMatrix<double> QK,
                     Eigen::SparseMatrix<double> Psi, Eigen::SparseMatrix<double> A,
                     int Ns, double tol, bool verbose = false) {
  // Bring in the spde matrices
  Eigen::SparseMatrix<double> Cmat     = Eigen::SparseMatrix<double> (spde["Cmat"]);
  Eigen::SparseMatrix<double> Gmat     = Eigen::SparseMatrix<double> (spde["Gmat"]);
  Eigen::SparseMatrix<double> GtCinvG     = Eigen::SparseMatrix<double> (spde["GtCinvG"]);
  int n_spde = Cmat.rows();
  int K = theta.size();
  K = (K - 1) / 2;
  int sig2_ind = 2*K;
  Eigen::SparseMatrix<double> AdivS2 = A / theta[sig2_ind];
  Eigen::SparseMatrix<double> Sig_inv = QK + AdivS2;
  SimplicialLLT<Eigen::SparseMatrix<double> > cholSigInv;
  cholSigInv.compute(Sig_inv);
  cholSigInv.analyzePattern(Sig_inv);
  if(verbose) {Rcout << "Initial theta: " << theta.transpose() << std::endl;}
  // Initialize everything
  Eigen::SparseMatrix<double> Xpsi = X * Psi;
  Eigen::SparseMatrix<double> Qk(n_spde, n_spde);
  Eigen::VectorXd XpsiY = Xpsi.transpose() * y;
  // Eigen::MatrixXd Avh = A * Vh;
  double yy = y.transpose() * y;
  // Regular fixed point updates
  // Eigen::VectorXd theta_new;
  // Eigen::VectorXd thetaDiff(2 * K + 1);
  // int Step = 0;
  // double eps = tol + 1;
  // // while(eps > tol | Step < 5) {
  // while(Step < 1) {
  //   theta_new = theta_fixpt(theta, A, QK, cholSigInv, XpsiY, Xpsi, Vh, Avh, y, yy, spde, tol);
  //   Step += 1;
  //   Rcout << "Step " << Step << " theta: " << theta_new.transpose() << std::endl;
  //   thetaDiff = theta_new - theta;
  //   // This uses the percent
  //   // thetaDiff = thetaDiff.array() / theta.array();
  //   // eps = thetaDiff.maxCoeff();
  //   // Rcout << "eps = " << eps << std::endl;
  //   // The following uses the sqrt of the squared norm to determine whether to stop
  //   eps = thetaDiff.squaredNorm();
  //   eps = std::sqrt(eps);
  //   theta = theta_new;
  // }
  // Using SQUAREM
  SquaremOutput SQ_result;
  SquaremDefault.tol = tol;
  SQ_result = theta_squarem2(theta, A, QK, cholSigInv, XpsiY, Xpsi, Ns, y, yy, spde, tol, verbose);
  theta= SQ_result.par;
  // Bring results together for output
  if(verbose) {Rcout << "Final theta: " << theta.transpose() << std::endl;}
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

