//
//  ztrti2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "trti2.h"

using LATL::TRTI2;

int ztrti2_(char &uplo,char &diag, int &n,complex<double> *A,int &ldA,int &info)
{
   info=TRTI2<double>(uplo,diag,n,A,ldA);
   if(info!=0)
   {
      info=-info;
      xerbla_("ZTRTI2 ",info);
   }
   return 0;
}
