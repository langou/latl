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

int dgeqrf_(int &m,int &n,complex<double> *A,int &ldA,complex<double> *tau,complex<double> *work,int &lwork,int &info)
{
   int nb=80;
   work[0]=n*nb;
   if(lwork<0)
      return 0;
   else
      nb=lwork/n;
   info=GEQRF<double>(m,n,A,ldA,tau,nb,work);
   return 0;
}
