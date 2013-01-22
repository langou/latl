//
//  dsymm.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "symm.h"

using latl::symm;

int dsymm_(char& side, char& uplo, int &m, int &n, double &alpha, double *A, int &ldA, double *B, int &ldB, double &beta, double *C, int &ldC)
{
   char str[]="DSYMM ";
   int info=-symm<double>(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   if(info>0)
      xerbla_(str,info);
   return 0;
}
