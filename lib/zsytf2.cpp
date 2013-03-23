//
//  zsytf2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/26/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "sytf2.h"
#include "latl.h"

using LATL::sytf2;
using LATL::int_t;

int zsytf2_(char &uplo,int &n,complex<double> *A,int &ldA,int *ipiv,int &info)
{
   int_t *IPIV=new int_t[n];
   bool *BSDV=new bool[n];
   info=sytf2<double>(uplo,n,A,ldA,IPIV,BSDV);
   for(int i=0;i<n;i++)
   {
      ipiv[i]=IPIV[i]+1;
      if(BSDV[i])
         IPIV[i]=-IPIV[i];
   }
   delete [] IPIV;
   delete [] BSDV;
   if(info<0)
      xerbla_("ZSYTF2  ",info);
   return 0;
}
