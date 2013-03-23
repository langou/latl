//
//  ctrmm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmm.h"

using LATL::trmm;

int ctrmm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB)
{
   int info=-trmm<float>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info!=0)
      xerbla_("CTRMM ",info);
   return 0;
}

