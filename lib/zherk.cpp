//
//  zherk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "herk.h"

using LATL::HERK;
using std::complex;

int zherk_(char& uplo, char& trans, int &n, int& k, double &alpha, complex<double> *A, int &ldA, double &beta, complex<double> *C, int &ldC)
{
   int info=-HERK<double>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("ZHERK ",info);
   return 0;
}
