//
//  cgeqrf.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "geqrf.h"

using LATL::GEQRF;

int cgeqrf_(int &m,int &n,complex<float> *A,int &ldA,complex<float> *tau,complex<float> *work,int &lwork,int &info)
{
   int nb=80;
   work[0]=n*nb;
   if(lwork<0)
      return 0;
   else
      nb=lwork/n;
   info=GEQRF<float>(m,n,A,ldA,tau,nb,work);
   return 0;
}
