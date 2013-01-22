//
//  srotmg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotmg.h"

using latl::rotmg;

int srotmg_(float& d1, float& d2, float& x1, float& y1, float *param)
{
   rotmg<float>(d1,d2,x1,y1,param);
   return 0;
}
