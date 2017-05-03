R:=Rscript
REPORT_DIR:=work/report
FIGURE_DIR= figures
CACHEDIR= cache
LMK=latexmk -pdf -f --interaction=nonstopmode -outdir=$(REPORT_DIR) -bibtex
.PHONY= all clean

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

artifacts/%.pdf: $(REPORT_DIR)/%.pdf
	cp $< $@

$(REPORT_DIR)/%.bib: doc/%.bib
	-mkdir -p $(REPORT_DIR)	
	cp $< $@

$(REPORT_DIR)/%.pdf: $(REPORT_DIR)/%.tex $(REPORT_DIR)/bibliography.bib
	$(LMK) $<

$(REPORT_DIR)/%.pdf: doc/%.tex $(REPORT_DIR)/bibliography.bib
	$(LMK) $<

include dependencies.makefile
