//
//  strsm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trsm.h"

using latl::trsm;

int strsm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, float& alpha, float *A, int &ldA, float *B, int &ldB)
{
   int info=-trsm<float>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info>0)
      xerbla_("STRSM ",info);
   return 0;
}
