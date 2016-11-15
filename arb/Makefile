INSTALL:=${HOME}/devel/install

GMP_LIB_DIR       = $(INSTALL)/lib
GMP_INCLUDE_DIR   = $(INSTALL)/include
MPFR_LIB_DIR      = $(INSTALL)/lib
MPFR_INCLUDE_DIR  = $(INSTALL)/include
FLINT_LIB_DIR     = $(INSTALL)/lib
FLINT_INCLUDE_DIR = $(INSTALL)/include/flint
ARB_LIB_DIR     = $(INSTALL)/lib
ARB_INCLUDE_DIR = $(INSTALL)/include/

CC=gcc
CFLAGS:=-O2 -ggdb -ansi -std=c99 -pedantic -Wall -funroll-loops
CFLAGS:=-ggdb -ansi -std=c99 -pedantic -Wall

LIBS=-L$(ARB_LIB_DIR) -L$(FLINT_LIB_DIR) -L$(GMP_LIB_DIR) -L$(MPFR_LIB_DIR) -larb -lflint -lmpfr -lgmp -lpthread -lm
INCS=-I$(CURDIR) -I$(GMP_INCLUDE_DIR) -I$(MPFR_INCLUDE_DIR) -I$(FLINT_INCLUDE_DIR) -I$(ARB_INCLUDE_DIR)
LINKLIBS=$(GMP_LIB_DIR)/libgmp.a $(MPFR_LIB_DIR)/libmpfr.a $(FLINT_LIB_DIR)/libflint.so -I$(ARB_INCLUDE_DIR)/libarb.a

LD_LIBRARY_PATH:=${CURDIR}:${FLINT_LIB_DIR}:${MPFR_LIB_DIR}:${GMP_LIB_DIR}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
export

BUILD:=_build
PROG   =
ifeq ($(shell uname -s),Darwin)
LIBRARY   = libajsec.dylib
else
LIBRARY   = libajsec.so
endif

SOURCES = $(wildcard *.c)

OBJS = $(patsubst %.c, $(BUILD)/%.o, $(SOURCES))
LOBJS = $(patsubst %.c, $(BUILD)/%.lo, $(SOURCES))
LIB_SOURCES = $(SOURCES) $(foreach dir, $(BUILD_DIRS), $(wildcard $(dir)/*.c))
LIB_OBJS = $(patsubst %.c, $(BUILD)/%.lo, $(LIB_SOURCES))

$(BUILD)/%.lo: %.c $(BUILD)
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

$(BUILD)/%.o: %.c $(BUILD)
	$(CC) -fPIC $(CFLAGS) $(INCS) -c $< -o $@

# programs
prog: $(LIBRARY)
	$(foreach prog, $(PROG), $(CC) $(CFLAGS) $(INCS) $(prog).c -o $(prog) $(LIBS);)

# library
library: $(LIBRARY)
$(LIBRARY): $(OBJS) $(LIB_OBJS)
	$(CC) -fPIC -shared $(LIB_OBJS) $(LIBS2) -o $(LIBRARY)

$(BUILD):
	mkdir -p $(BUILD)

