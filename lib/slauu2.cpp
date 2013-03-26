//
//  slauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauu2.h"

using LATL::LAUU2;

int slauu2_(char &uplo,int &n,float *A,int &ldA,int &info)
{
   info=LAUU2<float>(uplo,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("SLAUU2 ",info);
   }
   return 0;
}
