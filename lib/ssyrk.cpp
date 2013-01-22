//
//  ssyrk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using latl::syrk;

int ssyrk_(char& uplo, char& trans, int &n, int& k, float& alpha, float *A, int &ldA, float &beta, float *C, int &ldC)
{
   int info=-syrk<float>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("SSYRK ",info);
   return 0;
}
