//
//  strtri.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trtri.h"

using LATL::TRTRI;

int strtri_(char &uplo,char &diag, int &n,float *A,int &ldA,int &info)
{
   int nb=40;
   info=TRTRI<float>(uplo,diag,n,A,ldA,nb);
   if(info!=0)
   {
      info=-info;
      xerbla_("STRTRI ",info);
   }
   return 0;
}
