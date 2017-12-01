.PHONY= all clean
R:=Rscript
RNW_FILES=$(wildcard *.Rnw)
RENDER_FILES=$(addprefix render_,$(RNW_FILES:.Rnw=.tex))
TEX_FILES=$(wildcard *.tex) $(RENDER_FILES)
LMK=latexmk -pdf -xelatex --interaction=nonstopmode 

all: vignette.pdf
vignette.pdf: $(TEX_FILES)

clean: rnwclean

rnwclean: 
	-rm $(RENDER_FILES)


test:
	@echo $(TEX_FILES)

%.pdf: %.tex
	$(LMK) $<

render_%.tex: %.Rnw
	$R  -e "devtools::load_all()" \
		-e "knitr::opts_knit[['set']](root.dir = normalizePath('./'))" \
		-e "knitr::opts_chunk[['set']](echo=FALSE, warning=FALSE)" \
		-e "knitr::opts_chunk[['set']](cache.path='$(basename $(@F))/')" \
		-e "knitr::opts_chunk[['set']](fig.path='$(basename $(@F))/')" \
		-e "knitr::opts_chunk[['set']](fig.lp='fig:')" \
		-e "knitr::opts_chunk[['set']](results='asis')" \
		-e "knitr::opts_chunk[['set']](dpi=144, fig.width=4, fig.height=3)" \
		-e "knitr::knit('$<', output='$@')"
