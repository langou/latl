//
//  zhpr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hpr.h"

using LATL::HPR;
using std::complex;

int zhpr_(char& uplo, int &n, double &alpha, complex<double> *x, int& incx, complex<double> *A)
{
   int info=-HPR<double>(uplo,n,alpha,x,incx,A);
   if(info>0)
      xerbla_("ZHPR  ",info);
   return 0;
}
