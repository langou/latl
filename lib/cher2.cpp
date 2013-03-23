//
//  cher2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "her.h"

using LATL::her;
using std::complex;

int cher2_(char& uplo, int &n, complex<float> &alpha, complex<float> *x, int& incx, complex<float> *y, int& incy,complex<float> *A, int &ldA)
{
   int info=-her<float>(uplo,n,alpha,x,incx,y,incy,A,ldA);
   if(info>0)
      xerbla_("CHER2  ",info);
   return 0;
}
