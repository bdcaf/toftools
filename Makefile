
data/sample_ions.Rdata: dev_only/make_data.R
	Rscript '$<'

tmp/documented: $(wildcard *.R)
	R -e 'library(devtools);document()'
	touch $@

tofTools_0.2.tar.gz: $(wildcard *.R) tmp/documented roxygen
	R CMD build .


attributes:
	Rscript -e 'Rcpp::compileAttributes()'

roxygen: attributes
	Rscript -e 'roxygen2::roxygenise()'

build: tofTools_0.2.tar.gz data/sample_ions.Rdata
