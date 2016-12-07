TEX:= AbelJacobiMap.tex ComputationalAspects.tex ComputingAbelJacobi.tex DoubleExponentialIntegration.tex Intersections.tex Introduction.tex StrategyPeriodMatrix.tex SuperellipticCurves.tex
PDF:=$(TEX:.tex=.pdf)
%.pdf: %.tex macros.tex packages.tex
	pdflatex $<

default: $(PDF)

all: $(wildcard *.tex)
	pdflatex main.tex
