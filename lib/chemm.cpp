//
//  chemm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "hemm.h"

using std::complex;
using latl::hemm;

int chemm_(char& side, char& uplo, int &m, int &n, complex<float> &alpha, complex<float> *A, int &ldA, complex<float> *B, int &ldB, complex<float> &beta, complex<float> *C, int &ldC)
{
   int info=-hemm<float>(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("CHEMM ",info);
   return 0;
}
