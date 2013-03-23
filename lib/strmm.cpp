//
//  strmm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmm.h"

using LATL::TRMM;

int strmm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, float& alpha, float *A, int &ldA, float *B, int &ldB)
{
   int info=-TRMM<float>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info!=0)
      xerbla_("STRMM ",info);
   return 0;
}

