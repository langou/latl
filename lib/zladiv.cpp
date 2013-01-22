//
//  zladiv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "ladiv.h"

using latl::ladiv;

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__clang__)
complex<double> zladiv_(complex<double> &x, complex<double> &y)
{
   return ladiv<double>(x,y);
}
#else
int zladiv_(complex<double>& ret, complex<double> &x, complex<double> &y)
{
   ret=ladiv<double>(x,y);
   return 0;
}
#endif


