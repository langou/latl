//
//  dlauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/21/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using latl::lauum;

int dlauum_(char &uplo,int &n,double *A,int &ldA,int &info)
{
   int nb=80;
   info=lauum<double>(uplo,n,A,ldA,nb);
   if(info!=0)
   {
      info=-info;
      xerbla_("DLAUUM ",info);
   }
   return 0;
}
