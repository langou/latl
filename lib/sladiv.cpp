//
//  sladiv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "ladiv.h"

using LATL::ladiv;

int sladiv_(float &A, float &B, float &C, float &D, float &P, float &Q)
{
   ladiv<float>(A,B,C,D,P,Q);
   return 0;
}
