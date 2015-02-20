#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec loglog_lnlCpp(arma::vec theta, arma::mat y, arma::mat X) {
  
  int ncols = X.n_cols - 1;
  arma::vec beta = theta.subvec(0, ncols);
  int g = theta[ncols+1];
  arma::vec d = y.submat(0,0,y.n_rows-1,0); 
  arma::vec ti = y.submat(0,1,y.n_rows-1,1); 
  arma::vec ly = y.submat(0,2,y.n_rows-1,2); 
  arma::vec t0 = y.submat(0,3,y.n_rows-1,3); 

  arma::vec lambda = exp(-X * beta);
  int alpha = exp(-g);
  arma::vec lnFt = pow((lambda % ti), alpha);

  return(lnFt);
}


