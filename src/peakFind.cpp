#include <Rcpp.h>
using namespace Rcpp;

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
NumericVector readScales(NumericVector ys, NumericVector pos){
  //int n = ys.size();
  int k = pos.size();
  NumericVector out(k);

  out[0] = NA_REAL;
  //int sp = ceil(pos[0]);
  //int fin = ceil(pos[k-1]);

  for (int ip = 1; ip<k; ip++) {
	// is zero based in R so -1
	double left = pos[ip-1]-1;
	double right = pos[ip]-1;

	int li = ceil(left);
	int ri = floor(right);

	double count = 0;
	count += ys[floor(left)] * (  li - left);
	count -= ys[ri] * (1- (right -   ri));
	for (int is = li; is<=ri; is++) {
	  count += ys[is];
	}
	out[ip]=count;
  }
  return(out);
}




// [[Rcpp::export]]
NumericVector testRcp(int n){
  NumericVector out(n);

  for (int i=0; i<n; i++){
  	out[i] = i;
  }

  return(out);
}
