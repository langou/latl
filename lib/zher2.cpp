//
//  zher2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "her2.h"

using LATL::HER2;
using std::complex;

int zher2_(char& uplo, int &n, complex<double> &alpha, complex<double> *x, int& incx, complex<double> *y, int& incy,complex<double> *A, int &ldA)
{
   int info=-HER2<double>(uplo,n,alpha,x,incx,y,incy,A,ldA);
   if(info>0)
      xerbla_("ZHER2  ",info);
   return 0;
}
