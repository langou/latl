//
//  rotg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _rotg_h
#define _rotg_h

/// @file rotg.h Constructs Givens plane rotations.

#include <cmath>
#include "latl.h"

namespace latl 
{
   /// @brief Constructs real Givens plane rotations.
   ///
   /// For a point (x,y) in the plane, c and s are determined so that
   ///
   ///        [  c  s ] [ x ]  =  [ r ]
   ///        [ -s  c ] [ y ]     [ 0 ]
   /// where r = sqrt(x^2+y^2).
   /// @tparam real_t Floating point type.
   /// @param[in,out] x Real scalar, the x coordinate.
   /// @param[in,out] y Real scalar, the y coordinate.
   /// @param[out] c Real scalar, cos(theta).
   /// @param[out] s Real scalar, sin(theta).
   /// @ingroup ROT

   template <typename real_t>
   void rotg(real_t &x, real_t &y, real_t &c, real_t &s)
   {
      using std::abs;
      using std::sqrt;

      const real_t zero(0.0);
      const real_t one(1.0);
      real_t roe,r,scale,z,tx,ty;
      
      if(abs(x)>abs(y))
         roe=x;
      else
         roe=y;
      scale=abs(x)+abs(y);
      
      if(scale==zero)
      {
         c=one;
         s=zero;
         r=zero;
         z=zero;
      }
      else
      {
         tx=(x/scale);
         ty=(y/scale);
         if(roe>zero)
            r=scale*sqrt(tx*tx+ty*ty);
         else
            r=-scale*sqrt(tx*tx+ty*ty);
         c=x/r;
         s=y/r;
         if(abs(x)>abs(y))
            z=s;
         else if(c!=zero)
            z=one/c;
         else
            z=one;
      }
      x=r;
      y=z;
   }
   
   /// @brief Constructs complex Givens plane rotations.
   ///
   /// For a point (x,y) in the complex plane, c and s are determined so that
   ///
   ///        [  c        s ] [ x ]  =  [ r ]
   ///        [ -conj(s)  c ] [ y ]     [ 0 ]
   /// where r = (x/|x|) sqrt(|x|^2 + |y|^2) if |x|!=0, or r = y if |x|=0.
   /// @tparam real_t Floating point type.
   /// @param[in,out] x Complex scalar, the x coordinate.
   /// @param[in,out] y Complex scalar, the y coordinate.
   /// @param[out] c Real scalar, cos(theta).
   /// @param[out] s Complex scalar, sin(theta).
   /// @ingroup ROT

   template <typename real_t>
   void rotg(complex<real_t> &x,complex<real_t> y, real_t &c, complex<real_t> &s)
   {
      using std::conj;
      using std::sqrt;
      using std::abs;

      const real_t zero(0.0);
      const real_t one(1.0);
      complex<real_t> alpha;
      real_t norm,scale,tx,ty;
      if(abs(x)==zero)
      {
         c=zero;
         s=one;
         x=y;
      }
      else
      {
         scale=abs(x)+abs(y);
         tx=abs(x/scale);
         ty=abs(y/scale);
         norm=scale*sqrt(tx*tx+ty*ty);
         alpha=x/abs(x);
         c=abs(x)/norm;
         s=alpha*conj(y)/norm;
         x=alpha*norm;
      }
   }
}
#endif

