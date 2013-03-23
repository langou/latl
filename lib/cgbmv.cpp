//
//  cgbmv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gbmv.h"

using LATL::GBMV;
using std::complex;

int cgbmv_(char& trans, int &m, int &n, int& kL, int& kU, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *x, int& incx, complex<float> &beta, complex<float> *y, int& incy)
{
   int info=-GBMV<float>(trans,m,n,kL,kU,alpha,A,ldA,x,incx,beta,y,incy);
   if(info>0)
      xerbla_("CGBMV ",info);
   return 0;
}
