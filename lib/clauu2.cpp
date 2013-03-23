//
//  clauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using LATL::lauum;

int clauu2_(char &uplo,int &n,complex<float> *A,int &ldA,int &info)
{
   info=lauum<float>(uplo,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("CLAUU2 ",info);
   }
   return 0;
}
