//
//  sgeqrf.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "geqrf.h"

using LATL::GEQRF;

int sgeqrf_(int &m,int &n,float *A,int &ldA,float *tau,float *work,int &lwork,int &info)
{
   int nb=80;
   info=GEQRF<float>(m,n,A,ldA,tau,nb);
   return 0;
}
