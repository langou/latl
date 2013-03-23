//
//  cgerc.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "ger.h"

using LATL::gerc;
using std::complex;

int cgerc_(int &m, int &n, complex<float> &alpha, complex<float> *x, int& incx, complex<float> *y, int& incy, complex<float> *A, int &ldA)
{
   int info=-gerc<float>(m,n,alpha,x,incx,y,incy,A,ldA);
   if(info!=0)
      xerbla_("CGERC  ",info);
   return 0;
}
