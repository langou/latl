//
//  cherk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "herk.h"

using LATL::herk;
using std::complex;

int cherk_(char& uplo, char& trans, int &n, int& k, float& alpha, complex<float> *A, int &ldA, float &beta, complex<float> *C, int &ldC)
{
   int info=-herk<float>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("CHERK ",info);
   return 0;
}
