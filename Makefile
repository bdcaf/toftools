
tmp/documented: $(wildcard *.R)
	R -e 'library(devtools);document()'
	touch $@

tofTools_0.2.tar.gz: $(wildcard *.R) tmp/documented
	R CMD build .

build: tofTools_0.2.tar.gz
