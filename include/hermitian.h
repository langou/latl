//
//  hermitian.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _hermitian_h
#define _hermitian_h

/// @file hermitian.h Generates random complex Hermitian matrices.

#include <cstdlib>
#include <random>
#include "latl.h"

namespace LATL
{
   /// @brief Generates random Hermitian n-by-n complex matrix.
   /// @tparam real_t Floating point type.
   /// @return Pointer to matrix of type std::complex<real_t> containing n-by-n random matrix.
   /// @return nullptr if parameter invalid.
   /// @param n Number of columns of matrix.  n>0.
   /// @param s Seed for random number generator (optional).
   /// @ingroup MATGEN

   template<typename real_t> complex<real_t> *Hermitian(int_t n,uint32_t s=0)
   {
      using std::conj;
      if(n<1)
         return nullptr;

      complex<real_t> *M=new complex<real_t>[n*n];
      if(M)
      {
         std::random_device device;
         std::mt19937 generator(device());
         if(s>0)
            generator.seed(s);
         std::uniform_real_distribution<real_t> distribution(0,1);
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;i<j;i++)
            {
               M[i+j*n]=complex<real_t>(distribution(generator),distribution(generator));
               M[j+i*n]=conj(M[i+j*n]);
            }
            M[j+j*n]=distribution(generator);
         }
      }
      return M;
   }
}
#endif
