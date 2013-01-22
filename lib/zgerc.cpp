//
//  zgerc.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "ger.h"

using latl::gerc;
using std::complex;

int zgerc_(int &m, int &n, complex<double> &alpha, complex<double> *x, int& incx, complex<double> *y, int& incy, complex<double> *A, int &ldA)
{
   int info=-gerc<double>(m,n,alpha,x,incx,y,incy,A,ldA);
   if(info!=0)
      xerbla_("ZGERC  ",info);
   return 0;
}
