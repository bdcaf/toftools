data/sample_ions.Rdata: dev_only/make_data.R
	Rscript '$<'

tmp/documented: $(wildcard *.R) attributes
	R -e 'library(devtools);document()'
	touch $@

tofTools_0.2.tar.gz: $(wildcard *.R) tmp/documented 
	R CMD build .


attributes:
	Rscript -e 'Rcpp::compileAttributes()'

build: tofTools_0.2.tar.gz data/sample_ions.Rdata

install: data/sample_ions.Rdata attributes
	R CMD install --build .
