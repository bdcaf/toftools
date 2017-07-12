reports: $(addprefix artifacts/,report_warping.pdf)

$(REPORT_DIR)/%.pdf: doc/knitr_preamble.tex
$(REPORT_DIR)/vignette.pdf: doc/introduction.tex
$(REPORT_DIR)/report_warping.pdf: $(addprefix doc/,introduction.tex) \
	$(addprefix $(REPORT_DIR)/,warp_sample.tex)



