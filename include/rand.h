//
//  rand.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _rand_h
#define _rand_h

/// @file rand.h Generates random matrices.

#include <cstdlib>
#include <random>
#include "latl.h"

namespace LATL
{
   /// @brief Generates random m-by-n real matrix.
   /// @tparam real_t Floating point type.
   /// @return Pointer to matrix of type T containing m-by-n random matrix.
   /// @return NULL if parameter invalid.
   /// @param m Number of rows of matrix.  m>0.
   /// @param n Number of columns of matrix.  n>0.
   /// @ingroup MATGEN

   template<typename real_t> real_t *Rand(int_t m,int_t n)
   {
      if((n<1)||(m<1))
         return nullptr;

      real_t *M=new real_t[m*n];
      if(M)
      {
         std::random_device device;
         std::mt19937 generator(device());
         std::uniform_real_distribution<real_t> distribution(0,1);
         for(int_t i=0;i<m*n;i++)
            M[i]=distribution(generator);
      }
      return M;
   }

}
#endif
