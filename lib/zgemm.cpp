//
//  zgemm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "gemm.h"

using LATL::gemm;
using std::complex;

int zgemm_(char& transA, char& transB, int &m, int &n, int& k, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC)
{
   int info=-gemm<double>(transA,transB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("ZGEMM ",info);
   return 0;
}
