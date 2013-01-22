//
//  zsyrk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using std::complex;
using latl::syrk;

int zsyrk_(char& uplo, char& trans, int &n, int& k, complex<double> &alpha, complex<double> *A, int &ldA, complex<double> &beta, complex<double> *C, int &ldC)
{
   int info=-syrk<double>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("ZSYRK ",info);
   return 0;
}
