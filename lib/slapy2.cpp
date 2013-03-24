//
//  slapy2.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lapy2.h"

using LATL::LAPY2;

float slapy2_(float &x, float &y)
{
   return LAPY2<float>(x,y);
}


