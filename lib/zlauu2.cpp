//
//  zlauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using LATL::lauum;

int zlauu2_(char &uplo,int &n,complex<double> *A,int &ldA,int &info)
{
   info=lauum<double>(uplo,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("ZLAUU2 ",info);
   }
   return 0;
}
