# why
define VIGNETTE_OPTIONS
	%\VignetteIndexEntry{Acetone experiment warping}
	%\VignetteAuthor{Clemens Ager}
	%\VignetteKeyword{R}
	%\VignetteKeyword{package}
	%\VignetteKeyword{vignette}
	%\VignetteKeyword{example}
	%\VignetteKeyword{acetone}
	%\VignetteEngine{vignetteEngineMake::make}
	%\VignetteTangle{FALSE}
endef

THIS:=acetone_warp

$(THIS).pdf:
	cd $(THIS); make
	cd $(THIS); latexmk vignette.tex
	cp $(THIS)/vignette.pdf $@

include Makefile
