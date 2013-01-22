//
//  dtrmm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmm.h"

using latl::trmm;

int dtrmm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB)
{
   char str[]="DTRMM ";
   int info=-trmm<double>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info!=0)
      xerbla_(str,info);
   return 0;
}

