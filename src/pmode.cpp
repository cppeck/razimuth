#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//[[Rcpp::export]]
List pmode_cpp(arma::mat x, arma::mat ux){

  int In = x.n_rows;
  int Jn = ux.n_rows;
  arma::mat modeI(In, Jn, arma::fill::zeros);
  for (int i = 0; i < In; i++){
    for (int j = 0; j < Jn; j++){
      if(x(i,0) == ux(j,0) & x(i,1) == ux(j,1)) {
        modeI(i,j) = 1.0;
      } else {
        modeI(i,j) = 0.0;
      }
    }
  }

  return Rcpp::List::create(
    _["modeI"] = modeI);

}
