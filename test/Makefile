CC=clang++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL)

all: invert

float: clean
	@$(MAKE) all DEFS="-DREAL=float"

double: clean
	@$(MAKE) all DEFS="-DREAL=double"

ldouble: clean
	@$(MAKE) all DEFS="-DREAL=ldouble"

invert: invert.o
	$(CC) -o $@ invert.o

.cpp.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

clean:
	rm -f invert invert.o
