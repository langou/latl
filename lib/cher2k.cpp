//
//  cher2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "herk.h"

using LATL::HERK;
using std::complex;

int cher2k_(char& uplo, char& trans, int &n, int& k, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, float &beta, complex<float> *C, int &ldC)
{
   int info=-HERK<float>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("CHER2K ",info);
   return 0;
}
