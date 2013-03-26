//
//  dsyr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr2.h"

using LATL::SYR2;

int dsyr2_(char& uplo, int &n, double &alpha, double *x, int& incx, double *y, int& incy, double *A, int &ldA)
{
   int info=-SYR2<double>(uplo,n,alpha,x,incx,y,incy,A,ldA);
   if(info>0)
      xerbla_("DSYR2 ",info);
   return 0;
}
