//
//  zhemv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hemv.h"

using std::complex;
using LATL::hemv;

int zhemv_(char& uplo, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *x, int& incx, complex<double> &beta, complex<double> *y, int& incy)
{
   int info=-hemv<double>(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("ZHEMV ",info);
   return 0;
}
