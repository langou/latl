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

#if defined(FLOAT)
typedef float real_t;
#elif defined(DOUBLE)
typedef double real_t;
#elif defined(LDOUBLE)
typedef long double real_t;
#elif defined(MPREAL)
#include "mpreal.h"
using mpfr::mpreal;
typedef mpreal real_t;
#else
typedef double real_t;
#endif

int main(int argc, char** argv)
{
   using std::cout;
   using std::endl;
   using std::numeric_limits;
   using std::max;
#ifdef MPREAL
   int bits=53;
   if(argc>1)
      bits=max(2,atoi(argv[1]));
   mpreal::set_default_prec(bits);
#endif
   cout << endl;
   cout << "-------------------------------------" << endl;
   cout << "machine epsilon  = " << numeric_limits<real_t>::epsilon() << endl;
   cout << "minimum value    = " << numeric_limits<real_t>::min() << endl;
   cout << "maximum value    = " << numeric_limits<real_t>::max() << endl;
   cout << "-------------------------------------" << endl;
   cout << endl;
   return(0);
}
