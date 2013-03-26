//
//  dtrti2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trti2.h"

using LATL::TRTI2;

int dtrti2_(char &uplo,char &diag, int &n,double *A,int &ldA,int &info)
{
   info=TRTI2<double>(uplo,diag,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("DTRTI2 ",info);
   }
   return 0;
}
