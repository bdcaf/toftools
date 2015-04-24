library(Rcpp)

cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')

add(1,2,3)


cppFunction('double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}')

cppFunction('LogicalVector peakFind(NumericVector ys, int win, int min) {
  int n = ys.size();
  LogicalVector out(n);

  for(int i = win; i < n-win; ++i) {
    out[i] = ys[i] > min && ys[i]>ys[i-1] && ys[i]>ys[i+1];
  }
  return out;
}')
pf <- peakFind(sc, win=15, min=10)
plot(pf, type='l')
which(pf)
sc[pf]
