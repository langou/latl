//
//  genral.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _general_h
#define _general_h

/// @file general.h Generates random general matrices.

#include <cstdlib>
#include <random>
#include "latl.h"

namespace LATL
{
   /// @brief Generates random m-by-n real matrix.
   /// @tparam real_t Floating point type.
   /// @return Pointer to matrix of type real_t containing m-by-n random matrix.
   /// @return nullptr if parameter invalid.
   /// @param m Number of rows of matrix.  m>0.
   /// @param n Number of columns of matrix.  n>0.
   /// @param s Seed for random number generator (optional).
   /// @ingroup MATGEN

   template<typename real_t> real_t *General(int_t m,int_t n,uint32_t s=0)
   {
      if((n<1)||(m<1))
         return nullptr;

      real_t *M=new real_t[m*n];
      if(M)
      {
         std::random_device device;
         std::mt19937 generator(device());
         if(s>0)
            generator.seed(s);
         std::uniform_real_distribution<real_t> distribution(0,1);
         for(int_t i=0;i<m*n;i++)
            M[i]=distribution(generator);
      }
      return M;
   }

   /// @brief Generates random m-by-n complex matrix.
   /// @tparam real_t Floating point type.
   /// @return Pointer to matrix of type std::complex<real_t> containing m-by-n random matrix.
   /// @return nullptr if parameter invalid.
   /// @param m Number of rows of matrix.  m>0.
   /// @param n Number of columns of matrix.  n>0.
   /// @param s Seed for random number generator (optional).
   /// @ingroup MATGEN
   
   template<typename real_t> complex<real_t> *GeneralComplex(int_t m,int_t n,uint32_t s=0)
   {
      if((n<1)||(m<1))
         return nullptr;

      complex<real_t> *M=new complex<real_t>[m*n];
      if(M)
      {
         std::random_device device;
         std::mt19937 generator(device());
         if(s>0)
            generator.seed(s);
         std::uniform_real_distribution<real_t> distribution(0,1);
         for(int_t i=0;i<m*n;i++)
            M[i]=complex<real_t>(distribution(generator),distribution(generator));
      }
      return M;
   }
}
#endif
