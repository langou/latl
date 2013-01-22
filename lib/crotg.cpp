//
//  crotg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotg.h"

using latl::rotg;

int crotg_(complex<float>& a,complex<float>& b, float& c, complex<float>& s)
{
   rotg<float>(a,b,c,s);
   return 0;
}
