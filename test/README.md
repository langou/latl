This directory contains source for test codes using latl functions.  
The Makefile targets correspond to the precision to use in the test codes:

   float, double, ldouble (long double) use the built-in types

   real53, real128, real256, real512, real1024 use mpfr::real<> (see below)
    
   mpreal uses mpfr::mpreal (see below)
    
In order to use the multiprecision MPFR-based types, you will need the following:

   1. The GNU Multiprecision library (GMP) from http://gmplib.org
   2. The MPFR library from http://www.mpfr.org
   3. Pavel Holoborodko's mpfr::mpreal from http://www.holoborodko.com/pavel/mpfr
    
       - or -
       
   4. Christian Schneider's mpfr::real from http://chschneider.eu/programming/mpfr_real/
    
The mpfr::mpreal class is contained in the header file mpreal.h, and the precision
can be set at runtime.

The mpfr::real class is contained in the header file real.hpp, and the precision is
set at compile time as a template parameter (hence, the real53, real128, etc. targets).

The program lamch will print out the machine parameters for the selected precision.

The program invert will read a (text format) matrix from standard input and attempt
to compute its inverse, then print out the maximum relative residual out of
|| A*inv(A)-I || and || inv(A)*A-I ||.  When the mpreal precision is used, the
precision is set with the -P command line option.  Note that several types of matrices
can be inverted.  A general matrix is default, but symmetric, Hermitian, and positive
definite can also be selected (see invert --help), as well as complex.

Test matrices can be generated from the matgen program in ../matgen/.
    
