//
//  dsyr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr.h"

using LATL::SYR;

int dsyr_(char& uplo, int &n, double &alpha, double *x, int& incx, double *A, int &ldA)
{
   int info=-SYR<double>(uplo,n,alpha,x,incx,A,ldA);
   if(info>0)
      xerbla_("DSYR  ",info);
   return 0;
}
