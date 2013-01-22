//
//  larnv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/22/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _larnv_h
#define _larnv_h

/// @file larnv.h Returns a vector of random numbers from a uniform or normal distribution.

#include <algorithm>
#include <cmath>
#include <random>
#include "latl.h"

namespace latl
{
   /// @brief Returns a vector of n random real numbers from a uniform or normal distribution.
   ///
   /// Requires ISO C++ 2011 random number generators.
   /// @tparam real_t Floating point type.
   /// @param dist Specifies the distribution:
   ///
   ///        1: uniform (0,1)
   ///        2: uniform (-1,1)
   ///        3: normal (0,1)
   /// @param n Length of vector x.
   /// @param x Pointer to real vector of length n.
   /// @ingroup VEC
   
   template<typename real_t> 
   void larnv(int_t dist,int_t n,real_t *x)
   {
      std::random_device device;
      std::mt19937 generator(device());
      if(dist==1)
      {
         std::uniform_real_distribution<real_t> d1(0,1);
         for(int_t i=0;i<n;i++)
            x[i]=d1(generator);
      }
      else if(dist==2)
      {
         std::uniform_real_distribution<real_t> d2(-1,1);            
         for(int_t i=0;i<n;i++)
            x[i]=d2(generator);
      }
      else if(dist==3)
      {
         std::normal_distribution<real_t> d3(0,1);
         for(int_t i=0;i<n;i++)
            x[i]=d3(generator);
      }
   }

   /// @brief Returns a vector of n random complex numbers from a uniform or normal distribution.
   ///
   /// Requires ISO C++ 2011 random number generators.
   /// @tparam real_t Floating point type.
   /// @param dist Specifies the distribution of the random numbers:
   ///
   ///        1:  real and imaginary parts each uniform (0,1)
   ///        2:  real and imaginary parts each uniform (-1,1)
   ///        3:  real and imaginary parts each normal (0,1)
   ///        4:  uniformly distributed on the disc abs(z) < 1
   ///        5:  uniformly distributed on the circle abs(z) = 1
   /// @param n Length of vector x.
   /// @param x Pointer to complex vector of length n.
   /// @ingroup VEC
   
   template<typename real_t> 
   void larnv(int_t dist,int_t n,complex<real_t> *x)
   {
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t eight=8.0;
      const real_t twopi=eight*std::atan(one);
      std::random_device device;
      std::mt19937 generator(device());
      if(dist==1)
      {
         std::uniform_real_distribution<real_t> d1(0,1);
         for(int_t i=0;i<n;i++)
            x[i]=complex<real_t>(d1(generator),d1(generator));
      }
      else if(dist==2)
      {
         std::uniform_real_distribution<real_t> d2(-1,1);            
         for(int_t i=0;i<n;i++)
            x[i]=complex<real_t>(d2(generator),d2(generator));
      }
      else if(dist==3)
      {
         std::normal_distribution<real_t> d3(0,1);
         for(int_t i=0;i<n;i++)
            x[i]=complex<real_t>(d3(generator),d3(generator));
      }
      else if(dist==4)
      {
         std::uniform_real_distribution<real_t> d4(0,1);            
         for(int_t i=0;i<n;i++)
            x[i]=std::sqrt(d4(generator))*std::exp(complex<real_t>(zero,twopi*d4(generator)));
      }
      else if(dist==5)
      {
         std::uniform_real_distribution<real_t> d5(0,1);            
         for(int_t i=0;i<n;i++)
            x[i]=std::exp(complex<real_t>(zero,twopi*d5(generator)));
      }
   }
}

#endif
