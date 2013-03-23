//
//  zher.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "her.h"

using LATL::her;
using std::complex;

int zher_(char& uplo, int &n, double &alpha, complex<double> *x, int& incx, complex<double> *A, int &ldA)
{
   int info=-her<double>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("ZHER  ",info);
   return 0;
}
