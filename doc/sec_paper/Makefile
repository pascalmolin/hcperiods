TEX:= AbelJacobiMap.tex ComputationalAspects.tex ComputingAbelJacobi.tex NumericalIntegration.tex Intersections.tex Introduction.tex StrategyPeriodMatrix.tex SuperellipticCurves.tex ExamplesTimings.tex Outlook.tex
PDF:=$(TEX:.tex=.pdf)
%.pdf: %.tex macros.tex packages.tex
	pdflatex $<

default: $(PDF)

all: $(wildcard *.tex)
	pdflatex main.tex

IPE:= $(wildcard images/*.ipe)
images/%.pdf: images/%.ipe
	ipetoipe -pdf $< $@

ipe: $(IPE:.ipe=.pdf)

zip: $(TEX) main.tex macros.tex packages.tex main.bib $(wildcard images/*)
	zip article.zip $+
