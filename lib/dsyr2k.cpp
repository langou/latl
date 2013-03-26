//
//  dsyr2k.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "syr2k.h"

using LATL::SYR2K;

int dsyr2k_(char& uplo, char& trans, int &n, int& k, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC)
{
   int info=-SYR2K<double>(uplo,trans,n,k,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_("DSYR2K ",info);
   return 0;
}
