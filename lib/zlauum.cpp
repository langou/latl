//
//  zlauum.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/21/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lauum.h"

using LATL::LAUUM;

int zlauum_(char &uplo,int &n,complex<double> *A,int &ldA,int &info)
{
   int nb=80;
   info=LAUUM<double>(uplo,n,A,ldA,nb);
   if(info!=0)
   {
      info=-info;
      xerbla_("ZLAUUM ",info);
   }
   return 0;
}
