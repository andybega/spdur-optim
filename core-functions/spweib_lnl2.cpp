#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double spweib_lnlC(LogicalVector c, NumericVector pr0, NumericVector pr1,
  NumericVector st, NumericVector st0, NumericVector ln_ft) {
  int n = c.size();
  double lnl = 0;
  
  for(int i = 0; i < n; ++i) {
    if (c[i] == true) {
      lnl += log(pr0[i] + (pr1[i] * st[i])) - log(pr0[i] + (pr1[i] * st0[i]));
    } else {
      lnl += log(pr1[i]) + ln_ft[i]         - log(pr0[i] + (pr1[i] * st0[i]));
    }
  }
  
  return -lnl;
}


