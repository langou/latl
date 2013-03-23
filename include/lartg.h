//
//  lartg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/24/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lartg_h
#define _lartg_h

/// @file lartg.h Generates a plane rotation to eliminate second component.

#include <cmath>
#include <algorithm>
#include <limits>
#include "lapy2.h"
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
   ///        * if g=0, then c=1 and s=0.
   ///        * if f=0 and g!=0, then c=0 and s=1 without doing any computation
   /// If |f|>|g|, c>0.
   /// @tparam real_t Floating point type.
   /// @param f The first component to be rotated.
   /// @param g The second component to be rotated.
   /// @param c The cosine of the rotation [out].
   /// @param s the sine of the rotation [out].
   /// @param r The nonzero component of the rotated vector [out].
   /// @ingroup ROT
   
   template<typename real_t>
   void lartg(real_t f, real_t g, real_t &c, real_t &s, real_t &r)
   {
      using std::numeric_limits;
      using std::conj;
      using std::real;
      using std::imag;
      using std::pow;
      using std::log;
      using std::max;
      using std::sqrt;
      using std::abs;
      using std::trunc;
      using LATL::lapy2;
      
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t two=2.0;
      real_t safmin=numeric_limits<real_t>::min();
      real_t eps=numeric_limits<real_t>::epsilon();
      real_t base=numeric_limits<real_t>::radix;
      real_t p=trunc((log(safmin/eps)/log(base))/two);
      real_t SafeMin=pow(base,p);
      real_t SafeMax=one/SafeMin;
      if(g==zero)
      {
         c=one;
         s=zero;
         r=f;
      }
      else if(f==zero)
      {
         c=zero;
         s=one;
         r=g;
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
         if((abs(f)>abs(g))&&(c<zero))
         {
            c=-c;
            s=-s;
            r=-r;
         }
      }
   }

   /// @brief Generates a plane rotation of a complex vector (f,g) to eliminate g.
   /// generate a plane rotation so that
   ///
   ///        [       c  s ] [ f ]  =  [ r ]   
   ///        [ -conj(s) c ] [ g ]     [ 0 ]
   /// where c^2+|s|^2=1.
   /// This is a faster version of complex LATL::rotg with the following differences:
   ///
   ///        * f and g are unchanged on return.
   ///        * if g=0, then c=1 and s=0.
   ///        * if f=0, then c=0 and is chosen so that r is real.
   /// @tparam real_t Floating point type.
   /// @param f The first component to be rotated.
   /// @param g The second component to be rotated.
   /// @param c The cosine of the rotation [out].
   /// @param s The sine of the rotation [out].
   /// @param r The nonzero component of the rotated vector [out].
   /// @ingroup ROT
   
   template<typename real_t>
   void lartg(complex<real_t> f, complex<real_t> g, real_t &c, complex<real_t> &s, complex<real_t> &r)
   {
      using std::numeric_limits;
      using std::conj;
      using std::real;
      using std::imag;
      using std::pow;
      using std::log;
      using std::max;
      using std::sqrt;
      using std::abs;
      using std::trunc;
      using LATL::lapy2;

      const real_t zero(0.0);
      const real_t one(1.0);
      const real_t two(2.0);
      const complex<real_t> czero(0.0,0.0);
      real_t safmin=numeric_limits<real_t>::min();
      real_t eps=numeric_limits<real_t>::epsilon();
      real_t base=numeric_limits<real_t>::radix;
      real_t p=trunc((log(safmin/eps)/log(base))/two);
      real_t SafeMin=pow(base,p);
      real_t SafeMax=one/SafeMin;
      complex<real_t> F=f;
      complex<real_t> G=g;
      real_t scale=max(max(abs(real(F)),abs(imag(F))),max(abs(real(G)),abs(imag(G))));
      if((g==czero)&&(scale<=SafeMin))
      {
         c=one;
         s=czero;
         r=f;
      }
      else
      {
         int count=0;
         while(scale>=SafeMax)
         {
            count++;
            F*=SafeMin;
            G*=SafeMin;
            scale*=SafeMin;
         }
         while(scale<=SafeMin)
         {
            count--;
            F*=SafeMax;
            G*=SafeMax;
            scale*=SafeMax;
         }
         real_t F2=norm(F);
         real_t G2=norm(G);
         if(F2<=max(G2,one)*safmin)
         {
            if(f==czero)
            {
               c=zero;
               r=lapy2(real(g),imag(g));
               s=conj(G)/lapy2(real(G),imag(G));
            }
            else
            {
               real_t F2S=lapy2(real(F),imag(F));
               real_t G2S=sqrt(G2);
               c=F2S/G2S;
               complex<real_t> ff;
               if(max(abs(real(f)),abs(imag(f)))>one)
                  ff=f/lapy2(real(f),imag(f));
               else
                  ff=SafeMax*f/lapy2(real(f),imag(f));
               s=ff*conj(G)/G2S;
               r=c*f+s*g;
            }
         }
         else
         {
            real_t temp=sqrt(one+G2/F2);
            r=temp*F;
            c=one/temp;
            real_t d=F2+G2;
            s=conj(G)*(r/d);
            if(count>0)
               for(int i=0;i<count;i++)
                  r*=SafeMax;
            else if(count<0)
               for(int i=0;i<-count;i++)
                  r*=SafeMin;
         }
         
      }
   }
}

#endif
