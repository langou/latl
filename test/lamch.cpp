//
//  lamch.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <cmath>

#if defined(FLOAT)
typedef float real_t;
#elif defined(DOUBLE)
typedef double real_t;
#elif defined(LDOUBLE)
typedef long double real_t;
#elif defined(REAL)
#include "real.hpp"
typedef mpfr::real<REAL> real_t;
#else
typedef double real_t;
#endif

int main(int argc, char** argv)
{
   using std::cout;
   using std::endl;
   using std::numeric_limits;

   cout << endl;
   cout << "-------------------------------------" << endl;
   cout << "machine epsilon  = " << numeric_limits<real_t>::epsilon() << endl;
   cout << "minimum value    = " << numeric_limits<real_t>::min() << endl;
   cout << "maximum value    = " << numeric_limits<real_t>::max() << endl;
   cout << "-------------------------------------" << endl;
   cout << endl;
   
   return(0);
}
