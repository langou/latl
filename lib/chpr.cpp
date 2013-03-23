//
//  chpr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
#include "blas.h"
#include "hpr.h"

using LATL::HPR;
using std::complex;

int chpr_(char& uplo, int &n, float& alpha, complex<float> *x, int& incx, complex<float> *A)
{
   int info=-HPR<float>(uplo,n,alpha,x,incx,A);
   if(info>0)
      xerbla_("CHPR  ",info);
   return 0;
}
