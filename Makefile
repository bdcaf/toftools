R:=Rscript
REPORT_DIR:=work/report
FIGURE_DIR= figures
CACHEDIR= cache
LMK=latexmk -pdf -f --interaction=nonstopmode -outdir=$(REPORT_DIR)
.PHONY= all clean

.DEFAULT_GOAL := all

all: reports

-include dependencies.makefile

$(REPORT_DIR): 
	mkdir -p $@

$(REPORT_DIR)/%: | $(REPORT_DIR)

work/documented: $(wildcard *.R) attributes
	R -e 'library(devtools);document()'
	touch $@

tofTools_0.2.tar.gz: $(wildcard *.R) tmp/documented 
	R CMD build .

attributes:
	Rscript -e 'Rcpp::compileAttributes()'

build: tofTools_0.2.tar.gz data/sample_ions.Rdata

install: data/sample_ions.Rdata attributes
	R CMD install --build .

prep_libs: prep_base
prep_base:
	$R -e "install.packages(c('DEoptim','BB','MALDIquant','MASS','Matrix','PeakSegDP','baseline','data.table','dplyr','ggplot2','nloptr','ptw','purrr','tidyr','zoo'))"
prep_bioc:
	$R -e 'source("https://bioconductor.org/biocLite.R")' \
	   -e 'biocLite(c("rhdf5"), ask=F)'

artifacts/%.pdf: $(REPORT_DIR)/%.pdf
	cp $< $@

$(REPORT_DIR)/%.tex: doc/%.Rnw
	$R  -e "require(knitr)" \
		-e "knitr::opts_knit[['set']](root.dir = normalizePath('./'))" \
		-e "knitr::opts_chunk[['set']](cache.path='$(REPORT_DIR)/$(basename $(@F))/')" \
		-e "knitr::opts_chunk[['set']](fig.path='$(REPORT_DIR)/$(basename $(@F))/')" \
		-e "knitr::opts_chunk[['set']](fig.lp='fig:')" \
		-e "knitr::opts_chunk[['set']](echo=FALSE, warning=FALSE, results='asis', dpi=144, fig.width=4, fig.height=3)" \
		-e "knitr::knit('$<', output='$@')"

$(REPORT_DIR)/%.pdf: $(REPORT_DIR)/%.tex
	$(LMK) $<

$(REPORT_DIR)/%.tex: doc/%.tex
	$(LMK) $<

clean: work_clean
work_clean: 
	rm -rf work/*

almost_clean:
	find work -iname "*.tex" -delete
	find work -iname "*.aux" -delete
	find work -iname "*.bib" -delete

