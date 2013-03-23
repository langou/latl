//
//  zsyr2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using LATL::SYRK;
using std::complex;

int zsyr2k_(char& uplo, char& trans, int &n, int& k, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> *B, int &ldB, complex<double> &beta, complex<double> *C, int &ldC)
{
   int info=-SYRK<double>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("ZSYR2K ",info);
   return 0;
}
