//
//  srotg.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "blas.h"
#include "rotg.h"

using LATL::ROTG;

int srotg_(float& A, float& B, float& C, float& S)
{
   ROTG<float>(A,B,C,S);
   return 0;
}
