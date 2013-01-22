//
//  csyrk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using std::complex;
using latl::syrk;

int csyrk_(char& uplo, char& trans, int &n, int& k, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> &beta, complex<float> *C, int &ldC)
{
   int info=-syrk<float>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("CSYRK ",info);
   return 0;
}
