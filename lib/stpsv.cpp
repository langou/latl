//
//  stpsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tpsv.h"

using latl::tpsv;

int stpsv_(char& uplo, char& trans, char& diag, int &n, float *A, float *x, int& incx)
{
   int info=-tpsv<float>(uplo,trans,diag,n,A,x,incx);
   if(info>0)
      xerbla_("STPSV ",info);
   return 0;
}
