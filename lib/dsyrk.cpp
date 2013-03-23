//
//  dsyrk.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using LATL::SYRK;

int dsyrk_(char& uplo, char& trans, int &n, int& k, double &alpha, double *A, int &ldA, double &beta, double *C, int &ldC)
{
   int info=-SYRK<double>(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC);
   if(info>0)
      xerbla_("DSYRK ",info);
   return 0;
}
