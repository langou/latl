CC=clang++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL)

all: invert lamch

float: clean
	@$(MAKE) all DEFS="-DFLOAT"

double: clean
	@$(MAKE) all DEFS="-DDOUBLE"

ldouble: clean
	@$(MAKE) all DEFS="-DLDOUBLE"

mpreal: clean
	@$(MAKE) all DEFS="-DMPREAL" LIBS="-lgmp -lmpfr"

invert: invert.o
	$(CC) -o $@ invert.o $(LIBS)

lamch: lamch.o
	$(CC) -o $@ lamch.o $(LIBS)

.cpp.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

clean:
	rm -f invert invert.o lamch lamch.o
