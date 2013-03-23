//
//  drotg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotg.h"

using LATL::rotg;

int drotg_(double& A, double& B, double& C, double& S)
{
   rotg<double>(A,B,C,S);
   return 0;
}
