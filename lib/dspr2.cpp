//
//  dspr2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "spr2.h"

using LATL::SPR2;

int dspr2_(char& uplo, int &n, double &alpha, double *x, int& incx, double *y, int& incy, double *A)
{
   int info=-SPR2<double>(uplo,n,alpha,x,incx,y,incy,A);
   if(info>0)
      xerbla_("DSPR2 ",info);
   return 0;
}
