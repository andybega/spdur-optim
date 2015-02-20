#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double loglog_lnlCpp(arma::vec theta, arma::mat y, arma::mat X) {
  
  int ncols = X.n_cols - 1;
  arma::vec beta = theta.subvec(0, ncols);
  double g = theta[ncols+1];
  arma::vec d = y.submat(0,0,y.n_rows-1,0); 
  arma::vec ti = y.submat(0,1,y.n_rows-1,1); 
  arma::vec ly = y.submat(0,2,y.n_rows-1,2); 
  arma::vec t0 = y.submat(0,3,y.n_rows-1,3); 

  arma::vec lambda = exp(-X * beta);
  double alpha = exp(-g);
  arma::vec lnFt = log(alpha) + alpha*log(lambda) + (alpha-1)*log(ti) - 2 * log(1 + pow(lambda % ti, alpha));
  arma::vec lnSt = -log(1 + pow(lambda % ti, alpha));
  arma::vec lnSt0 = -log(1 + pow(lambda % t0, alpha));

  arma::vec cens(lambda.size());
  for(int i=0; i<cens.size(); ++i) {
  	if ( d[i]==1 & ly[i]==0 | d[i]==0 ) {
  		cens[i] = lnSt[i] - lnSt0[i];
  	} else {
  		cens[i] = 0;
  	}
  }

  arma::vec nocens(lambda.size());
  for(int i=0; i<nocens.size(); ++i) {
  	if ( d[i]==1 & ly[i]==1 ) {
  		nocens[i] = lnFt[i] - lnSt0[i];
  	} else {
  		nocens[i] = 0;
  	}
  }

  double logL = sum(cens + nocens);
  return(logL);
}


