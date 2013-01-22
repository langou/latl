//
//  zher2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "herk.h"

using latl::herk;
using std::complex;

int zher2k_(char& uplo, char& trans, int &n, int& k, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, double &beta, complex<double> *C, int &ldC)
{
   int info=-herk<double>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("ZHER2K ",info);
   return 0;
}
