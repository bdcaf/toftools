read: $(REPORT_DIR)/report_warping.pdf
	open -a Skim.app  $<
reports: $(addprefix artifacts/,report_warping.pdf)

$(REPORT_DIR)/%.pdf: doc/knitr_preamble.tex doc/acronyms.tex
$(REPORT_DIR)/vignette.pdf: doc/introduction.tex
$(REPORT_DIR)/report_warping.pdf: $(addprefix doc/,introduction.tex) \
	$(addprefix $(REPORT_DIR)/,warp_sample.tex warp_insample.tex)



