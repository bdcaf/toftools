#include <Rcpp.h>
using namespace Rcpp;

//' find peaks within a window
//' 
//' Not really adjusted by my needs
// [[Rcpp::export]]
LogicalVector peakFind(NumericVector ys, int win, int min) {
  int n = ys.size();
  LogicalVector out(n);

  for(int i = win; i < n-win; ++i) {
	out[i] = ys[i] > min && ys[i]>ys[i-1] && ys[i]>ys[i+1];
  }
  return out;
}




// [[Rcpp::export]]
NumericVector testRcp(int n){
  NumericVector out(n);

  for (int i=0; i<n; i++){
  	out[i] = i;
  }

  return(out);
}
