read: $(REPORT_DIR)/report_warping.pdf
	open -a Skim.app  $<

reports: $(addprefix artifacts/,report_warping.pdf case_study_short_breaths.pdf)

$(REPORT_DIR)/%.pdf: $(addprefix doc/,knitr_preamble.tex acronyms.tex bibliography.bib) 

$(REPORT_DIR)/report_warping.pdf: $(addprefix doc/,introduction.tex) \
	$(addprefix $(REPORT_DIR)/,warp_sample.tex warp_insample.tex)

$(REPORT_DIR)/case_study_short_breaths.pdf: $(wildcard doc/case60_*.tex) \
  $(patsubst doc/%.Rnw,$(REPORT_DIR)/%.tex,$(wildcard doc/case60_*.Rnw))

data/sample_ions.Rdata: dev_only/make_data.R
	Rscript '$<'
