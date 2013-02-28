CC=clang++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL)

all: invert lamch

float:
	@$(MAKE) all DEFS="-DFLOAT"

double:
	@$(MAKE) all DEFS="-DDOUBLE"

ldouble:
	@$(MAKE) all DEFS="-DLDOUBLE"

real53:
	@$(MAKE) all DEFS="-DREAL=53" LIBS="-lgmp -lmpfr"

real128:
	@$(MAKE) all DEFS="-DREAL=128" LIBS="-lgmp -lmpfr"

real256:
	@$(MAKE) all DEFS="-DREAL=256" LIBS="-lgmp -lmpfr"

real512:
	@$(MAKE) all DEFS="-DREAL=512" LIBS="-lgmp -lmpfr"

real1024:
	@$(MAKE) all DEFS="-DREAL=1024" LIBS="-lgmp -lmpfr"

invert: invert.o
	$(CC) -o $@ invert.o $(LIBS)

lamch: lamch.o
	$(CC) -o $@ lamch.o $(LIBS)

.cpp.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

clean:
	rm -f invert invert.o lamch lamch.o
