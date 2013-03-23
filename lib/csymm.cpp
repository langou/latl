//
//  csymm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
///

#include "blas.h"
#include "symm.h"

using LATL::SYMM;
using std::complex;

int csymm_(char& side, char& uplo, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC)
{
   int info=-SYMM<float>(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("CSYMM ",info);
   return 0;
}
