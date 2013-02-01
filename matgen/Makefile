####
#### Builds matgen matrix generator.  NOTE: requires C++11 random number generator (gcc 4.7.2 works)
####

CC=/usr/local/bin/g++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL)

all: matgen

float: clean
	@$(MAKE) all DEFS="-DREAL=float"

double: clean
	@$(MAKE) all DEFS="-DREAL=double"

ldouble: clean
	@$(MAKE) all DEFS="-DREAL=ldouble"

matgen: matgen.o
	$(CC) -o $@ matgen.o

.cpp.o:
	$(CC) -c $(CFLAGS) $(DEFS) $<

clean:
	rm -f matgen matgen.o
