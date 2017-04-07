OUT:=$(wildcard out/*.out)

callgraphs: $(subst out,pdf,$(OUT))

pdf/%.pdf: out/%.out | pdf
	gprof2dot -f callgrind $< | dot -Tpdf -o $@

pdf:
	mkdir -p pdf

sync:
	rsync --progress -a epsilon:git/hcperiods/arb/profile/out .
