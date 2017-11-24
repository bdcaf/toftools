# why
define VIGNETTE_OPTIONS
  %\VignetteAuthor{Clemens Ager}
  %\VignetteKeyword{R}
  %\VignetteKeyword{package}
  %\VignetteKeyword{vignette}
  %\VignetteKeyword{Example}
  %\VignetteEngine{vignetteEngineMake::make}
  %\VignetteTangle{FALSE}
endef

THIS:=acetone_warp

$(THIS).pdf:
	cd $(THIS); latexmk vignette.tex
	cp $(THIS)/vignette.pdf $@

include Makefile
