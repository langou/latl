//
//  cgerc.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "geru.h"

using LATL::GERU;
using std::complex;

int cgeru_(int &m, int &n, complex<float> &alpha, complex<float> *x, int& incx, complex<float> *y, int& incy, complex<float> *A, int &ldA)
{
   int info=-GERU<float>(m,n,alpha,x,incx,y,incy,A,ldA);
   if(info!=0)
      xerbla_("CGERU  ",info);
   return 0;
}
