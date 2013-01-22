//
//  xerbla.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include <iostream>
#include <cstdlib>
using namespace std;

int xerbla_(const char *name, int &info)
{
   cout << " ** On entry to " << name;
   cout << " parameter number " << info;
   cout << " had an illegal value." << endl;
   exit(EXIT_FAILURE);
   return 0;
}
