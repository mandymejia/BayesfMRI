#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a number by two
//'
//' @param x a single integer
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}
