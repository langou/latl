//
//  stbsv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "tbsv.h"

using LATL::TBSV;

int stbsv_(char& uplo, char& trans, char& diag, int &n, int& k, float *A, int &ldA, float *x, int& incx)
{
   int info=-TBSV<float>(uplo,trans,diag,n,k,A,ldA,x,incx);
   if(info>0)
      xerbla_("STBSV ",info);
   return 0;
}
