//
//  ctrtri.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trtri.h"

using LATL::TRTRI;

int ctrtri_(char &uplo,char &diag, int &n,complex<float> *A,int &ldA,int &info)
{
   info=TRTRI<float>(uplo,diag,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("CTRTRI ",info);
   }
   return 0;
}
