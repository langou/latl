//
//  strti2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trti2.h"

using latl::trti2;

int strti2_(char &uplo,char &diag, int &n,float *A,int &ldA,int &info)
{
   info=trti2<float>(uplo,diag,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("STRTI2 ",info);
   }
   return 0;
}
