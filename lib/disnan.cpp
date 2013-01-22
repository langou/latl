//
//  disnan.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include <cmath>

using std::isnan;

int disnan_(double &din)
{
   return isnan(din);
}
