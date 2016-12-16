WHO=$(shell whoami)
ifeq ($(WHO), molin)
PREFIX=${HOME}/install
else ifeq ($(WHO), neurohr)
PREFIX=${HOME}/neurohr/superelliptic
else
PREFIX=/usr/local
endif
EXEC=export LD_LIBRARY_PATH=$(CURDIR):${PREFIX}/lib;
EXEEXT=

ifeq ($(shell uname -s),Darwin)
LIBRARY   = libajsec.dylib
else
LIBRARY   = libajsec.so
endif

INCS=-I$(CURDIR) -I/usr/include -I${PREFIX}/include -I${PREFIX}/include/flint
LIBS=-L$(CURDIR) -L/usr/lib -L${PREFIX}/lib -lajsec -larb -lflint -lmpfr -lgmp -lm -lpthread
LIBS2=-L$(CURDIR) -L/usr/lib -L${PREFIX}/lib -larb -lflint -lmpfr -lgmp -lm -lpthread

CC=gcc
CXX=g++
AR=ar
LDCONFIG=true

CFLAGS=-fno-common -ansi -pedantic -Wall -O2 -funroll-loops -g -mpopcnt -std=c99
ABI_FLAG=
PIC_FLAG=-fPIC
EXTRA_SHARED_FLAGS=
#-install_name /usr/local/lib/libarb-1.1.1.dylib -compatibility_version 1.1 -current_version 1.1.1

QUIET_CC  = @echo '   ' CC  ' ' $@;
QUIET_CXX = @echo '   ' CXX ' ' $@;
QUIET_AR  = @echo '   ' AR  ' ' $@;

AT=@

export

HEADERS = $(wildcard *.h)

SOURCES = $(wildcard *.c)
OBJS = $(patsubst %.c, build/%.o, $(SOURCES))
LOBJS = $(patsubst %.c, build/%.lo, $(SOURCES))

EXMP_SOURCES = $(wildcard examples/*.c)
EXMPS = $(patsubst %.c, %, $(EXMP_SOURCES))

TEST_SOURCES = $(wildcard test/*.c)
TESTS = $(patsubst %.c, build/%$(EXEEXT), $(TEST_SOURCES))

all: library

quiet: library

library: $(LIBRARY)

test: library $(TESTS)
	export LD_LIBRARY_PATH=${PREFIX}/lib
	$(AT)$(foreach prog, $(TESTS), $(EXEC) $(prog) || exit $$?;)

verbose:
	$(MAKE) AT= QUIET_CC= QUIET_CXX= QUIET_AR=

clean:
	rm -f $(OBJS) $(LOBJS) $(TESTS) $(PROFS) $(EXMPS) $(wildcard $(LIBRARY)*)
	rm -rf build

build:
	mkdir -p build

build/test:
	mkdir -p build/test

build/%.lo: %.c $(HEADERS) | build
	$(QUIET_CC) $(CC) $(PIC_FLAG) $(CFLAGS) $(INCS) -c $< -o $@;

build/%.o: %.c $(HEADERS) | build
	$(QUIET_CC) $(CC) $(CFLAGS) $(INCS) -c $< -o $@;

build/test/%$(EXEEXT): test/%.c $(HEADERS) $(LIBRARY) | build/test
	$(QUIET_CC) $(CC) $(CFLAGS) $(INCS) $< -o $@ $(LIBS)

$(LIBRARY): $(LOBJS) $(HEADERS) | build
	$(CC) $(ABI_FLAG) -shared $(EXTRA_SHARED_FLAGS) $(LOBJS) -o $(LIBRARY) $(LDFLAGS) $(LIBS2);

.PHONY: library clean examples check test distclean all valgrind
