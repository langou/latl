//
//  lartgs.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lartgs_h
#define _lartgs_h

/// @file lartgs.h Generates a plane rotation for a QR iteration.

#include <cmath>
#include "lamch.h"
#include "lartgp.h"
#include "latl.h"

namespace LATL
{
   /// @brief Generates a plane rotation for a QR iteration.
   ///
   /// Generates a plane rotation designed to introduce a bulge in
   /// Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD
   /// problem. If x and y are the top-row entries, and sigma is the shift,
   /// the computed c and s define a plane rotation satisfying
   ///
   ///        [  c  s ] [ x^2 - sigma ]  =  [ r ],
   ///        [ -s  c ] [    x * y    ]     [ 0 ]
   /// with r nonnegative.  If x^2-sigma and x*y are 0, then the
   /// rotation is by pi/2.
   /// @tparam real_t Floating point type.
   /// @param x The first diagonal element of an upper bidiagonal matrix.
   /// @param y The first super-diagonal element of an upper bidiagonal matrix.
   /// @param sigma The shift.
   /// @param c The cosine of the rotation.
   /// @param s The sine of the rotation.
   /// @ingroup COMP
   
   template<typename real_t>
   void LARTGS(real_t x, real_t y, real_t sigma, real_t &c, real_t &s)
   {
      using LATL::LARTGP;
      using std::abs;
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t thresh=LAMCH<real_t>('E');
      real_t W,Z,R,S;
      
      if(((sigma==zero)&&(abs(x)<thresh))||((abs(x)==sigma)&&(y==zero)))
      {
         Z=zero;
         W=zero;
      }
      else if(sigma==zero)
      {
         if(x>zero)
         {
            Z=-x;
            W=-y;
         }
         else
         {
            Z=-x;
            W=-y;
         }
      }
      else if(abs(x)<thresh)
      {
         Z=-sigma*sigma;
         W=zero;
      }
      else
      {
         S=(x>zero)?one:-one;
         Z=S*(abs(x)-sigma)*(S+sigma/x);
         W=s*y;
      }
      LARTGP<real_t>(W,Z,s,c,R);
   }
}

#endif
