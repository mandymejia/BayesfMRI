#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;
using namespace std;

/*
 B of size n1 x n2
 Set A(i:i+n1,j:j+n2) = B (update)
 */
// void setSparseBlock_update(SparseMatrix<double,0,int>* A,int i, int j, SparseMatrix<double,0,int>& B)
// {
//   for (int k=0; k<B.outerSize(); ++k) {
//     for (SparseMatrix<double,0,int>::InnerIterator it(B,k); it; ++it) {
//       A->coeffRef(it.row()+i, it.col()+j) = it.value();
//     }
//   }
// }

//' Get the prewhitening matrix using the AR coefficients and variance
//'
//' @param avg_AR a matrix with dimensions V by p, where V is the number of data
//'   locations and p is the AR order
//' @param nTime (integer) the length of the time series that is being prewhitened
//' @param avg_var a vector of length V containing the residual variances of the
//'   AR model
Eigen::SparseMatrix<double> getPWmatrix(Eigen::MatrixXd avg_AR, int nTime, Eigen::VectorXd avg_var) {
  int nVox = avg_AR.rows();
  int p = avg_AR.cols(), diffIJ;
  int NT = nVox * nTime;
  double sqrt_var, sqrt_prec;
  Eigen::SparseMatrix<double> PWmatrix(NT,NT), spSqrtInv(nTime,nTime);
  Eigen::SparseMatrix<double> spInv_v(nTime,nTime);
  Eigen::internal::BandMatrix<double> halfInv_v(nTime,nTime,p,p);
  Eigen::MatrixXd final_Inv_v(nTime,nTime), sqrtInv(nTime,nTime), eDinv_v(nTime,nTime);
  Eigen::MatrixXd Inv_v(nTime,nTime), hInv_v(nTime,nTime), thInv_v(nTime,nTime), Dinv_v(nTime,nTime);
  // Eigen::MatrixXd lowerInv_v(nTime, nTime), upperInv_v(nTime,nTime);
  // Eigen::DiagonalMatrix<double, Eigen::Dynamic> Dinv_v(nTime), eDinv(nTime);
  nVox = 1;
  for(int v=0; v<nVox; v++) {
    sqrt_var = sqrt(avg_var(v));
    sqrt_prec = 1/sqrt_var;
    //// Using BandMatrix
    Dinv_v.diagonal().setConstant(sqrt_prec);
    halfInv_v.diagonal().setConstant(1);
    for(int k=1; k<=p; k++) {
      // Inv_v.diagonal(k).setConstant(-avg_AR(v,k-1)); //superdiagonals
      halfInv_v.diagonal(-k).setConstant(-avg_AR(v,k-1)); //subdiagonals
    }
    hInv_v = halfInv_v.toDenseMatrix();
    thInv_v = hInv_v.transpose();
    Inv_v =  hInv_v * thInv_v;
    // Rcout << Inv_v.block(0,0,10,10) << std::endl;
    //// End using BandMatrix
    //// Without using BandMatrix
    // Now fill the off-diagonals with the negative coefficient values
    // for(int i; i<nTime; i++) {
    //   lowerInv_v(i,i) = 1;
    //   Dinv_v.diagonal()[i] = sqrt_prec;
    //   for(int j; j<i; j++) {
    //     diffIJ = i - j;
    //     for(int k; k<p; k++) {
    //       if(diffIJ == k) {
    //         lowerInv_v(i,j) = -avg_AR(v,k);
    //       }
    //     }
    //   }
    // }
    // upperInv_v = lowerInv_v.transpose();
    // Inv_v = lowerInv_v * upperInv_v;
    //// End not using BandMatrix
    final_Inv_v = Dinv_v * Inv_v * Dinv_v;
    // spInv_v = final_Inv_v.sparseView();
    // Rcout << spInv_v.block(0,0,10,10) << std::endl;
    // Eigen::EigenSolver<SparseMatrix<double>> ei;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ei(final_Inv_v);
    // if (ei.info() != Eigen::Success) {
      // Rcout << "Negative eigenvalues!" << std::endl;
    // }
    Eigen::VectorXd d2inv = ei.eigenvalues().real();
    Eigen::MatrixXd eVec = ei.eigenvectors().real();
    // std::vector<std::tuple<double, Eigen::VectorXd>> ei_vec_and_vals;
    // for(int i=0;i<d2inv.size();i++) {
    //   std::tuple<double,Eigen::VectorXd> vec_and_val(d2inv(i),eVec.col(i));
    //   ei_vec_and_vals.push_back(vec_and_val);
    // }
    // std::sort(ei_vec_and_vals.begin(),ei_vec_and_vals.end(),
    //           [&](const std::tuple<double, Eigen::VectorXd>& a, const std::tuple<double, Eigen::VectorXd>& b) -> bool{
    //             return std::get<0>(a) < std::get<0>(b);
    //           });
    // Rcout << "Unsorted eigenvalues: " << d2inv.transpose().reverse() << std::endl;
    // Rcout << "The first three eigenvectors are " << eVec.block(0,0,nTime,3) << std::endl;
    // Rcout << "The first eigenvector is " << eVec.transpose() << std::endl;
    // Rcout << d2inv << std::endl;
    // Eigen::VectorXd dinv = d2inv.sqrt();
    // Rcout << "d2inv.array().pow(.5) :" << d2inv.array().pow(.5) << std::endl;
    eDinv_v.diagonal() = d2inv.reverse();
    // for(int i;i<nTime;i++){
    //   eDinv_v(i,i) = sqrt(d2inv.reverse()(i));
    // }
    // Rcout << eDinv_v.block(0,0,10,10) << std::endl;
    // sqrtInv = eVec * eDinv_v.array().sqrt() * eVec.transpose();
    // Rcout << sqrtInv.block(0,0,10,10) << std::endl;
    // spSqrtInv = sqrtInv.sparseView();
    // setSparseBlock_update(&PWmatrix, v*nTime, v*nTime, spSqrtInv);
  }
  return PWmatrix;
}
