//
// GEMM timing test: to use vecLib on Mac OS, compile with
//
// /usr/local/bin/c++ -DDOUBLE -D__latl_cblas -std=c++11 -I../include -I /System/Library/Frameworks/vecLib.framework/Headers -lblas blas.cpp timer.cpp -o blas
//
//  or
// /usr/local/bin/c++ -DFLOAT -D__latl_cblas -std=c++11 -I../include -I /System/Library/Frameworks/vecLib.framework/Headers -lblas blas.cpp timer.cpp -o blas

#include <iostream>
#include <gemm.h>
#include <load.h>
#include <print.h>
#include <cstdlib>
#include "timer.h"

#if defined(FLOAT)
typedef float Real;
#elif defined(DOUBLE)
typedef double Real;
#elif defined(LDOUBLE)
typedef long double Real;
#elif defined(REAL)
#include "real.hpp"
typedef mpfr::real<REAL> Real;
#elif defined(MPREAL)
#include "mpreal.h"
typedef mpfr::mpreal Real;
#else
typedef double Real;
#endif
typedef complex<Real> Complex;


using namespace std;
using namespace LATL;

double ran()
{
   return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}

int main(int argc, char **argv)
{
   long n=1000;

   if(argc>1)
      n=atoi(argv[1]);

   Real *A=new Real[n*n];
   Real *B=new Real[n*n];
   Real *C=new Real[n*n];


   for(int i=0;i<n*n;i++)
   {
      A[i]=ran();
      B[i]=ran();
      C[i]=ran();
   }

   Timer T1;
   T1.Start();
   GEMM<Real>('N','N',n,n,n,1.0,A,n,B,n,1.0,C,n);
   T1.Stop();

   double flops=2.0*n*n*n;
   cout << flops/T1.Time() << endl;

}