//
//  sgemm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemm.h"

using LATL::GEMM;

int sgemm_(char& transA, char& transB, int &m, int &n, int& k, float& alpha, float *A, int &ldA, float *B, int &ldB, float &beta, float *C, int &ldC)
{
   int info=-GEMM<float>(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("SGEMM ",info);
   return 0;
}
