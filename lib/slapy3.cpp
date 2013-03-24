//
//  slapy3.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "lapy3.h"

using LATL::LAPY3;

float slapy3_(float &x, float &y, float &z)
{
   return LAPY3<float>(x,y,z);
}


