//
//  dtrsm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trsm.h"

using LATL::trsm;

int dtrsm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB)
{
   int info=-trsm<double>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info>0)
      xerbla_("DTRSM ",info);
   return 0;
}
