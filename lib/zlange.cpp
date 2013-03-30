//
//  zlange.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lange.h"

using LATL::LANGE;

double zlange_(char &norm,int &m,int &n,complex<double> *A,int &ldA,double *work)
{
   return LANGE<double>(norm,m,n,A,ldA);
}
