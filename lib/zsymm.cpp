//
//  zsymm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "symm.h"

using LATL::SYMM;
using std::complex;

int zsymm_(char& side, char& uplo, int &m, int &n, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC)
{
   int info=-SYMM<double>(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("ZSYMM ",info);
   return 0;
}
