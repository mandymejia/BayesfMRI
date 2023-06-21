#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

//' Get the prewhitening matrix for a single data location
//'
//' @param AR_coeffs a length-p vector where p is the AR order
//' @param nTime (integer) the length of the time series that is being prewhitened
//' @param avg_var a scalar value of the residual variances of the AR model
//' 
// [[Rcpp::export(.getSqrtInvCpp)]]
Eigen::SparseMatrix<double> getSqrtInvCpp(Eigen::VectorXd AR_coeffs, int nTime, double avg_var) {
  double sqrt_var = sqrt(avg_var);
  int p = AR_coeffs.size();
  double sqrt_prec = 1/sqrt_var;
  Eigen::VectorXd Dinv_v(nTime);
  for(int i=0; i<nTime;i++){Dinv_v(i) = sqrt_prec;}
  Eigen::MatrixXd matDinv_v = Dinv_v.asDiagonal();
  Eigen::MatrixXd halfInv_v(nTime,nTime);
  halfInv_v = MatrixXd::Zero(nTime,nTime);
  halfInv_v.diagonal().setConstant(1);
  for(int j=0;j<nTime;j++){
    for(int k=1;k<=p;k++) {
      if(j+k >nTime - 1){break;}
      halfInv_v(j+k,j) = -1 * AR_coeffs(k-1);
    }
  }
  Eigen::MatrixXd Inv_v = halfInv_v * halfInv_v.transpose();
  Eigen::MatrixXd final_Inv_v = matDinv_v * Inv_v * matDinv_v;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ei(final_Inv_v);
  Eigen::VectorXd d2inv = ei.eigenvalues().real();
  Eigen::MatrixXd eVec = ei.eigenvectors().real();
  Eigen::MatrixXd eDinv_v(nTime,nTime);
  eDinv_v = MatrixXd::Zero(nTime,nTime);
  for(int i=0;i<nTime;i++){
    eDinv_v(i,i) = sqrt(d2inv.reverse()(i));
  }
  Eigen::MatrixXd revEvec = eVec.rowwise().reverse();
  Eigen::MatrixXd sqrtInv = revEvec * eDinv_v * revEvec.transpose();
  Eigen::MatrixXd out(nTime,nTime);
  out = Eigen::MatrixXd::Zero(nTime,nTime);
  for(int j=0; j<nTime;j++){
    for(int k=-1*(p+1);k<=p+1;k++) {
      if(j + k < 0){continue;}
      if(j + k > nTime - 1){continue;}
      out(j+k,j) = sqrtInv(j+k,j);
    }
  }
  Eigen::SparseMatrix<double> final_out = out.sparseView(1e-8,1);
  return final_out;
}
