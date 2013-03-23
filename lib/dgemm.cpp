//
//  dgemm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemm.h"

using LATL::GEMM;

int dgemm_(char& transA, char& transB, int &m, int &n, int& k, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC)
{
   int info=-GEMM<double>(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("DGEMM ",info);
   return 0;
}
