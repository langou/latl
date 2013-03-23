//
//  dsymv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "symv.h"

using LATL::symv;

int dsymv_(char& uplo, int &n, double &alpha, double *A, int &ldA, double *x, int& incx, double &beta, double *y, int& incy)
{
   int info=-symv<double>(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
   if(info!=0)
      xerbla_("DSYMV ",info);
   return 0;
}
