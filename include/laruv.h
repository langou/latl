//
//  laruv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _laruv_h
#define _laruv_h

/// @file laruv.h Returns a vector of uniform random real numbers.

#include <random>
#include "latl.h"

namespace latl
{
   /// @brief Returns a vector of uniform random real numbers on (0,1).
   ///
   /// Requires the ISO C++ 2011 standard random number generator.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param x Pointer to real vector of length n.
   /// @ingroup VEC
   
   template<typename real_t>
   void laruv(int_t n, real_t *x)
   {
      std::random_device device;
      std::mt19937 generator(device());
      std::uniform_real_distribution<real_t> distribution(0,1);
      for(int_t i=0;i<n;i++)
         x[i]=distribution(generator);
   }
}
#endif
