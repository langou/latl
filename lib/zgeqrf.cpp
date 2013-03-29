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

int zgeqrf_(int &m,int &n,complex<double> *A,int &ldA,complex<double> *tau,complex<double> *work,int &lwork,int &info)
{
   int nb=80;
   info=GEQRF<double>(m,n,A,ldA,tau,nb);
   return 0;
}
