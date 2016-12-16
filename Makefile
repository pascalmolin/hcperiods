all: pdf arb

pdf:
	$(MAKE) -C doc/sec_paper all

arb:
	$(MAKE) -C arb

test:
	$(MAKE) -C arb test
