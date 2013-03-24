//
//  slauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using LATL::LAUUM;

int slauu2_(char &uplo,int &n,float *A,int &ldA,int &info)
{
   info=LAUUM<float>(uplo,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("SLAUU2 ",info);
   }
   return 0;
}
