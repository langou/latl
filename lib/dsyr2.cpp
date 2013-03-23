//
//  dsyr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr.h"

using LATL::syr;

int dsyr2_(char& uplo, int &n, double &alpha, double *x, int& incx, double *y, int& incy, double *A, int &ldA)
{
   int info=-syr<double>(uplo,n,alpha,x,incx,y,incy,A,ldA);
   if(info>0)
      xerbla_("DSYR2 ",info);
   return 0;
}
