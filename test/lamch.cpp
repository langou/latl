//
//  lamch.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <lamch.h>

#if defined(FLOAT)
typedef float real_t;
#elif defined(DOUBLE)
typedef double real_t;
#elif defined(LDOUBLE)
typedef long double real_t;
#elif defined(REAL)
#include "real.hpp"
typedef mpfr::real<REAL> real_t;
#elif defined(MPREAL)
#include "mpreal.h"
typedef mpfr::mpreal real_t;
#else
typedef double real_t;
#endif

int main(int argc, char** argv)
{
   using std::cout;
   using std::endl;
   using LATL::LAMCH;
#ifdef MPREAL
   int prec=53;
   using std::atoi;
   if(argc>1)
      prec=atoi(argv[1]);
   mpfr::mpreal::set_default_prec(prec);
#endif
   cout << "LAMCH('E')=" << LAMCH<real_t>('E') << endl;
   cout << "LAMCH('P')=" << LAMCH<real_t>('P') << endl;
   cout << "LAMCH('S')=" << LAMCH<real_t>('S') << endl;
   cout << "LAMCH('U')=" << LAMCH<real_t>('U') << endl;
   cout << "LAMCH('O')=" << LAMCH<real_t>('O') << endl;
   return(0);
}
