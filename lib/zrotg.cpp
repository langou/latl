//
//  zrotg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotg.h"

using LATL::ROTG;

int zrotg_(complex<double>& a,complex<double>& b, double& c, complex<double>& s)
{
   ROTG<double>(a,b,c,s);
   return 0;
}
