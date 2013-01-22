//
//  dlapy3.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lapy3.h"

using latl::lapy3;

double dlapy3_(double &x, double &y, double &z)
{
   return lapy3<double>(x,y,z);
}


