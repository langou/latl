//
//  dladiv.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include "ladiv.h"

using LATL::LADIV;

int dladiv_(double &A, double &B, double &C, double &D, double &P, double &Q)
{
   LADIV<double>(A,B,C,D,P,Q);
   return 0;
}


