CC=clang++
CXX=clang++
LATL=../include
CXXFLAGS=$(DEFS) -I$(LATL) -O3 -Wall -pedantic -std=c++11

TESTSUITE = 
PROGRAMS = invert lamch gebd2

all: $(PROGRAMS) $(TESTSUITE)

float: clean
	@$(MAKE) all DEFS="-DFLOAT" LDLIBS="-lblas"

double: clean
	@$(MAKE) all DEFS="-DDOUBLE" LDLIBS="-lblas"

ldouble: clean
	@$(MAKE) all DEFS="-DLDOUBLE"

real53: clean
	@$(MAKE) all DEFS="-DREAL=53" LDLIBS="-lgmp -lmpfr"

real128: clean
	@$(MAKE) all DEFS="-DREAL=128" LDLIBS="-lgmp -lmpfr"

real256: clean
	@$(MAKE) all DEFS="-DREAL=256" LDLIBS="-lgmp -lmpfr"

real512: clean
	@$(MAKE) all DEFS="-DREAL=512" LDLIBS="-lgmp -lmpfr"

real1024: clean
	@$(MAKE) all DEFS="-DREAL=1024" LDLIBS="-lgmp -lmpfr"

mpreal: clean
	@$(MAKE) all DEFS="-DMPREAL" LDLIBS="-lgmp -lmpfr"

gebd2: gebd2.o timer.o
	$(CXX) $^ -o $@

clean:
	rm -f *.o $(PROGRAMS) $(TESTSUITE)
