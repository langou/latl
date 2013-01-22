//
//  dlauu2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauu2.h"

using latl::lauu2;

int dlauu2_(char &uplo,int &n,double *A,int &ldA,int &info)
{
   info=lauu2<double>(uplo,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("DLAUU2 ",info);
   }
   return 0;
}
