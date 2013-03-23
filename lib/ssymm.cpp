//
//  ssymm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "symm.h"

using LATL::symm;

int ssymm_(char& side, char& uplo, int &m, int &n, float& alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC)
{
   int info=-symm<float>(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("SSYMM ",info);
   return 0;
}
