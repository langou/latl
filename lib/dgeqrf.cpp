//
//  dgeqrf.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "geqrf.h"

using LATL::GEQRF;

int dgeqrf_(int &m,int &n,double *A,int &ldA,double *tau,double *work,int &lwork,int &info)
{
   int nb=80;
   info=GEQRF<double>(m,n,A,ldA,tau,nb);
   return 0;
}
