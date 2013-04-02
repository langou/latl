CC=/usr/local/bin/c++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL) -I/System/Library/Frameworks/vecLib.framework/Headers

all: invert lamch blas gmres

float: clean
	@$(MAKE) all DEFS="-DFLOAT -D__latl_cblas" LIBS="-lblas"

double: clean
	@$(MAKE) all DEFS="-DDOUBLE -D__latl_cblas" LIBS="-lblas"

ldouble: clean
	@$(MAKE) all DEFS="-DLDOUBLE"

real53: clean
	@$(MAKE) all DEFS="-DREAL=53" LIBS="-lgmp -lmpfr"

real128: clean
	@$(MAKE) all DEFS="-DREAL=128" LIBS="-lgmp -lmpfr"

real256: clean
	@$(MAKE) all DEFS="-DREAL=256" LIBS="-lgmp -lmpfr"

real512: clean
	@$(MAKE) all DEFS="-DREAL=512" LIBS="-lgmp -lmpfr"

real1024: clean
	@$(MAKE) all DEFS="-DREAL=1024" LIBS="-lgmp -lmpfr"

mpreal: clean
	@$(MAKE) all DEFS="-DMPREAL" LIBS="-lgmp -lmpfr"

invert: invert.o
	$(CC) -o $@ invert.o $(LIBS)

lamch: lamch.o
	$(CC) -o $@ lamch.o $(LIBS)

blas: blas.o timer.o
	$(CC) -o $@ blas.o timer.o $(LIBS) $(CBLAS)

gmres: gmres.o
	$(CC) -o $@ gmres.o $(LIBS)

.cpp.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

clean:
	rm -f *.o invert lamch blas gmres
