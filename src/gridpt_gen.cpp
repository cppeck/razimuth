// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec priorUpts(arma::mat gridpts, arma::mat obspt, arma::vec radius) {
  
  int M = gridpts.n_rows;
  int J = obspt.n_rows;
  arma::mat dist(M,J);
  arma::mat keep(M,J);
  arma::vec save(M);
  
  arma::mat distTmp(M,J);
  for (int j = 0; j < J; j++){
    for (int m = 0; m < M; m++) {
      distTmp(m,j) = sum(arma::pow((obspt.row(j) - gridpts.row(m)),2.0));
    }
  }
  dist = arma::sqrt(distTmp);
  for (int j = 0; j < J; j++){
    for (int m = 0; m < M; m++) {
      if(dist(m,j) <= radius[j]){
        keep(m,j) = 1.0;
      } else {
        keep(m,j) = 0.0;
      }
    }
  }
  for (int j = 0; j < J; j++){
    for (int m = 0; m < M; m++){
      if(sum(keep.row(m)) > 0.0){
        save[m] = 1;
      } else {
        save[m] = 0;
      }
    }
  }
  
  return save;
}