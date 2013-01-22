####
#### Makefile to build legacy BLAS+LAPACK interface library
####

CC=clang++
LATL=../include
CFLAGS=-O3 -Wall -pedantic -std=c++11 -I$(LATL)
LIBTOOL=ar cr
LIB=liblatl.a

####
#### BLAS objects
####

OBJB=\
saxpy.o isamax.o sasum.o scopy.o sdot.o snrm2.o srot.o srotg.o sscal.o sswap.o sdsdot.o srotm.o srotmg.o \
daxpy.o idamax.o dasum.o dcopy.o ddot.o dnrm2.o drot.o drotg.o dscal.o dswap.o dsdot.o drotm.o drotmg.o \
caxpy.o icamax.o scasum.o ccopy.o cdotc.o cdotu.o scnrm2.o csscal.o cscal.o csrot.o crotg.o cswap.o \
zaxpy.o izamax.o dzasum.o zcopy.o zdotc.o zdotu.o dznrm2.o zdscal.o zscal.o zdrot.o zrotg.o zswap.o \
sgemv.o sgbmv.o ssymv.o ssbmv.o sspmv.o strmv.o stbmv.o stpmv.o strsv.o stbsv.o stpsv.o sger.o ssyr.o sspr.o ssyr2.o sspr2.o \
dgemv.o dgbmv.o dsymv.o dsbmv.o dspmv.o dtrmv.o dtbmv.o dtpmv.o dtrsv.o dtbsv.o dtpsv.o dger.o dsyr.o dspr.o dsyr2.o dspr2.o \
cgemv.o cgbmv.o chemv.o chbmv.o chpmv.o ctrmv.o ctbmv.o ctpmv.o ctrsv.o ctbsv.o ctpsv.o cgerc.o cgeru.o cher.o chpr.o cher2.o chpr2.o \
zgemv.o zgbmv.o zhemv.o zhbmv.o zhpmv.o ztrmv.o ztbmv.o ztpmv.o ztrsv.o ztbsv.o ztpsv.o zgerc.o zgeru.o zher.o zhpr.o zher2.o zhpr2.o \
sgemm.o ssymm.o ssyrk.o ssyr2k.o strmm.o strsm.o \
dgemm.o dsymm.o dsyrk.o dsyr2k.o dtrmm.o dtrsm.o \
cgemm.o csymm.o csyrk.o csyr2k.o ctrmm.o ctrsm.o chemm.o cherk.o cher2k.o \
zgemm.o zsymm.o zsyrk.o zsyr2k.o ztrmm.o ztrsm.o zhemm.o zherk.o zher2k.o

####
#### LAPACK objects
####

OBJL=xerbla.o dlamch.o slamch.o disnan.o sisnan.o slange.o dlange.o clange.o zlange.o sladiv.o dladiv.o cladiv.o zladiv.o slapy2.o dlapy2.o slapy3.o dlapy3.o slauu2.o dlauu2.o clauu2.o zlauu2.o slauum.o dlauum.o clauum.o zlauum.o

####
#### rules and targets
####

all: 
	@$(MAKE) lib OBJ="$(OBJB) $(OBJL)"

blas: 
	@$(MAKE) lib OBJ="$(OBJB)"

lib: $(LIB)

$(OBJB): blas.h

$(OBJL): lapack.h

$(LIB): $(OBJ)
	$(LIBTOOL) $@ $(OBJ)

.cpp.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f $(OBJB) $(OBJL) $(LIB)
