//
//  clange.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lange.h"

using latl::lange;

float clange_(char &norm,int &m,int &n,complex<float> *A,int &ldA,float *work)
{
   return lange<float>(norm,m,n,A,ldA,work);
}
