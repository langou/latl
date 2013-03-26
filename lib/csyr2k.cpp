//
//  csyr2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
#include "blas.h"
#include "syr2k.h"

using LATL::SYR2K;
using std::complex;

int csyr2k_(char& uplo, char& trans, int &n, int& k, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC)
{
   int info=-SYR2K<float>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("CSYR2K ",info);
   return 0;
}
