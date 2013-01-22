//
//  ztrmm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "trmm.h"

using latl::trmm;

int ztrmm_(char& side, char& uplo, char& trans, char& diag, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB)
{
   int info=-trmm<double>(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB);
   if(info!=0)
      xerbla_("ZTRMM ",info);
   return 0;
}

