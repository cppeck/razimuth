#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//[[Rcpp::export]]
arma::vec vhf_cc_mu(arma::vec cp, arma::vec thetalist, arma::mat thetatmp, int distid) {

	int M = thetatmp.n_rows;
    int v = cp.n_elem;

	// Initialize storage
	arma::vec num(M);
	int Z = thetalist.n_elem;

	// Loop
	if (distid == 1){
	  if (v == 1){
	    for (int m = 0; m < M; m++) {
	      double tmpSum = 0;
	      for (int z = 0; z < Z; z++){
	        tmpSum += cp[0] * cos(thetalist[z] - thetatmp(m,z));
	      }
	      num[m] = tmpSum;
	    }
	  } else {
	    for (int m = 0; m < M; m++) {
	      double tmpSum = 0;
	      for (int z = 0; z < Z; z++){
	        tmpSum += cp[z] * cos(thetalist[z] - thetatmp(m,z));
	      }
	      num[m] = tmpSum;
	    }
	  }
	}

	if (distid == 2){
	  if (v == 1){
	    for (int m = 0; m < M; m++) {
	      double tmpSum = 0;
	      for (int z = 0; z < Z; z++){
	        tmpSum += -log(1.0 + pow(cp[0], 2.0) - 2.0 * cp[0] * cos(thetalist[z] - thetatmp(m,z)));
	      }
	      num[m] = tmpSum;
	    }
	  } else {
	    for (int m = 0; m < M; m++) {
	      double tmpSum = 0;
	      for (int z = 0; z < Z; z++){
	        tmpSum += -log(1.0 + pow(cp[z], 2.0) - 2.0 * cp[z] * cos(thetalist[z] - thetatmp(m,z)));
	      }
	      num[m] = tmpSum;
	    }
	  }
	}

	double C = num.max();
	arma::vec adjTmp = num - C;
	arma::vec logProb = adjTmp - log(sum(exp(adjTmp)));
	arma::vec prob = exp(logProb);

	return prob;

}
