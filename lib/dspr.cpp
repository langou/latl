//
//  dspr.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "spr.h"

using LATL::spr;

int dspr_(char& uplo, int &n, double &alpha, double *x, int& incx, double *A)
{
   int info=-spr<double>(uplo,n,alpha,x,incx,A);
   if(info>0)
      xerbla_("DSPR  ",info);
   return 0;
}
