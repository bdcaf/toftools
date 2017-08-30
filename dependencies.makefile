read: $(REPORT_DIR)/report_warping.pdf
	open -a Skim.app  $<

reports: $(addprefix artifacts/,report_warping.pdf)

REPORT_DEP=$(addprefix doc/,knitr_preamble.tex acronyms.tex bibliography.bib) 
$(REPORT_DIR)/vignette.pdf: doc/introduction.tex
$(REPORT_DIR)/report_warping.pdf: $(addprefix doc/,introduction.tex) \
	$(addprefix $(REPORT_DIR)/,warp_sample.tex warp_insample.tex) $(REPORT_DEP)

data/sample_ions.Rdata: dev_only/make_data.R
	Rscript '$<'
