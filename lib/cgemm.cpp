//
//  cgemm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemm.h"

using latl::gemm;
using std::complex;

int cgemm_(char& transA, char& transB, int &m, int &n, int& k, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC)
{
   int info=-gemm<float>(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("CGEMM ",info);
   return 0;
}
