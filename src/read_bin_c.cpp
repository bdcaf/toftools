#include <Rcpp.h>
using namespace Rcpp;



//' calculates the number of counts between indices
//' 
//' @param ys spectrum vector
//' @param pos position vector
//' @export
//' @return vector of length pos with the summed counts (first is 0)
// [[Rcpp::export]]
NumericVector read_bin_c(NumericVector ys, NumericVector pos){
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
