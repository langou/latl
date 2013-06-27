This directory contains source code to build a LAPACK+BLAS compatible static
library callable from FORTRAN.  The BLAS portion is complete, but only a few
LAPACK functions are currently implemented.  

The library has been tested using GNU g++/gfortran 4.7.2 and higher and 
Intel icc/ifort 13.0 and higher on Linux and MacOSX.  Note that mixing the
C++ and FORTRAN compilers may break the FORTRAN compatibility, mainly due
to the convention used for complex-values functions.

The C++ header file blas.h contains prototypes for the FORTRAN BLAS functions,
and lapack.h contains (some of the) prototypes for the FORTRAN LAPACK functions.
Note that these may also be used to call the FORTRAN code from C++.
