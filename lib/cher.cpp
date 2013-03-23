//
//  cher.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "her.h"

using LATL::HER;
using std::complex;

int cher_(char& uplo, int &n, float& alpha, complex<float> *x, int& incx, complex<float> *A, int &ldA)
{
   int info=-HER<float>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("CHER  ",info);
   return 0;
}
