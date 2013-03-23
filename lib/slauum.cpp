//
//  slauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/21/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using LATL::lauum;

int slauum_(char &uplo,int &n,float *A,int &ldA,int &info)
{
   int nb=80;
   info=lauum<float>(uplo,n,A,ldA,nb);
   if(info!=0)
   {
      info=-info;
      xerbla_("SLAUUM ",info);
   }
   return 0;
}
