//
//  ssyr2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr2k.h"

using LATL::SYR2K;

int ssyr2k_(char& uplo, char& trans, int &n, int& k, float& alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC)
{
   int info=-SYR2K<float>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("SSYR2K ",info);
   return 0;
}
