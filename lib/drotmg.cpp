//
//  drotmg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotmg.h"

using latl::rotmg;

int drotmg_(double& d1, double& d2, double& x1, double& y1, double *param)
{
   rotmg<double>(d1,d2,x1,y1,param);
   return 0;
}
