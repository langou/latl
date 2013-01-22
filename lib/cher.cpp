//
//  cher.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "her.h"

using latl::her;
using std::complex;

int cher_(char& uplo, int &n, float& alpha, complex<float> *x, int& incx, complex<float> *A, int &ldA)
{
   int info=-her<float>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("CHER  ",info);
   return 0;
}
