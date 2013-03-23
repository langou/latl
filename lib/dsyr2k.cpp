//
//  dsyr2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syrk.h"

using LATL::syrk;

int dsyr2k_(char& uplo, char& trans, int &n, int& k, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC)
{
   int info=-syrk<double>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("DSYR2K ",info);
   return 0;
}
