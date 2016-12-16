PREFIX=/usr/local
EXEEXT=

ifeq ($(shell uname -s),Darwin)
LIBRARY   = libajsec.dylib
else
LIBRARY   = libajsec.so
endif

BASEDIR=/home/neurohr/superelliptic

#INCS=-I$(CURDIR) -I/usr/local/include -I/usr/local/include/flint -I
INCS=-I$(CURDIR) -I$(BASEDIR)/include
LIBS=-L$(CURDIR) -L$(BASEDIR)/lib -lajsec -larb -lflint -lmpfr -lgmp -lm -lpthread
LIBS2=-L$(CURDIR) -L$(BASEDIR)/lib -larb -lflint -lmpfr -lgmp -lm -lpthread

CC=gcc
CXX=g++
AR=ar
LDCONFIG=true

CFLAGS=-fno-common -ansi -pedantic -Wall -O2 -funroll-loops -g -mpopcnt
ABI_FLAG=
PIC_FLAG=-fPIC
EXTRA_SHARED_FLAGS=
#-install_name /usr/local/lib/libarb-1.1.1.dylib -compatibility_version 1.1 -current_version 1.1.1

DLPATH=DYLD_LIBRARY_PATH
DLPATH_ADD=$(CURDIR):/usr/local/lib
EXTENSIONS=
EXTRA_BUILD_DIRS=

ifdef $(DLPATH)
  $(DLPATH):=$($(DLPATH)):$(DLPATH_ADD)
else
  $(DLPATH):=$(DLPATH_ADD)
endif

LIBDIR=lib

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
	$(AT)$(foreach prog, $(TESTS), $(prog) || exit $$?;)

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

build/test/%$(EXEEXT): test/%.c $(HEADERS) | build/test
	$(QUIET_CC) $(CC) $(CFLAGS) $(INCS) $< -o $@ $(LIBS)

$(LIBRARY): $(LOBJS) $(HEADERS) | build
	$(CC) $(ABI_FLAG) -shared $(EXTRA_SHARED_FLAGS) $(LOBJS) -o $(LIBRARY) $(LDFLAGS) $(LIBS2);

.PHONY: library clean examples check test distclean all valgrind
