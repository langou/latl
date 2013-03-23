//
//  lartgp.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lartgp_h
#define _lartgp_h

/// @file lartgp.h Generates a plane rotation of a real vector.

#include <cmath>
#include <algorithm>
#include <limits>
#include "latl.h"

namespace LATL
{
   /// @brief Generates a plane rotation of a real vector (f,g) to eliminate g.
   /// generate a plane rotation so that
   ///
   ///        [  c  s ] [ f ]  =  [ r ]   
   ///        [ -s  c ] [ g ]     [ 0 ]
   /// where c^2+s^2=1.
   /// This is a slower, more accurate version of LATL::rotg with the following other differences:
   ///
   ///        * f and g are unchanged on return.
   ///        * if g=0, then c=1 or c=-1 and s=0.
   ///        * if f=0 and g!=0, then c=0 and s=1 or s=-1.
   ///        * the sign is chosen so that r is nonnegative.
   /// @tparam real_t Floating point type.
   /// @param f The first component to be rotated.
   /// @param g The second component to be rotated.
   /// @param c The cosine of the rotation.
   /// @param s the sine of the rotation.
   /// @param r The nonzero component of the rotated vector.
   /// @ingroup ROT
   
   template<typename real_t>
   void lartgp(real_t f, real_t g, real_t &c, real_t &s, real_t &r)
   {
      using std::abs;
      using std::sqrt;
      using std::log;
      using std::max;
      using std::numeric_limits;
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t two=2.0;
      real_t safmin=numeric_limits<real_t>::min();
      real_t eps=numeric_limits<real_t>::epsilon();
      real_t base=numeric_limits<real_t>::radix;
      real_t p=(log(safmin/eps)/log(base))/two;
      real_t SafeMin=pow(base,p);
      real_t SafeMax=one/SafeMin;
      if(g==zero)
      {
         c=(f>=zero)?one:-one;
         s=zero;
         r=abs(f);
      }
      else if(f==zero)
      {
         c=zero;
         s=(g>=zero)?one:-one;
         r=abs(g);
      }
      else
      {
         real_t F=f;
         real_t G=g;
         real_t scale=max(abs(F),abs(G));
         int count=0;
         if(scale>=SafeMax)
         {
            while(scale>=SafeMax)
            {
               count++;
               F*=SafeMin;
               G*=SafeMin;
               scale=max(abs(F),abs(G));
            }
            r=sqrt(F*F+G*G);
            c=F/r;
            s=G/r;
            for(int i=0;i<count;i++)
               r*=SafeMax;
         }
         else if(scale<=SafeMin)
         {
            while(scale<=SafeMin)
            {
               count++;
               F*=SafeMax;
               G*=SafeMax;
               scale=max(abs(F),abs(G));
            }
            r=sqrt(F*F+G*G);
            c=F/r;
            s=G/r;
            for(int i=0;i<count;i++)
               r*=SafeMin;
         }
         else
         {
            r=sqrt(F*F+G*G);
            c=F/r;
            s=G/r;
         }
         if(r<zero)
         {
            c=-c;
            s=-s;
            r=-r;
         }
      }
   }
}
   
#endif
