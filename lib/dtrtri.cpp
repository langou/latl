//
//  dtrtri.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trtri.h"

using LATL::TRTRI;

int dtrtri_(char &uplo,char &diag, int &n,double *A,int &ldA,int &info)
{
   info=TRTRI<double>(uplo,diag,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("DTRTRI ",info);
   }
   return 0;
}
