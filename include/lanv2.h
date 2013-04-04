//
//  lanv2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanv2_h
#define _lanv2_h

/// @file lanv2.h Computes the Schur factorization of a real 2-by-2 matrix.

#include <cmath>
#include <algorithm>
#include "lapy2.h"
#include "lamch.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes the Schur factorization of a real 2-by-2 matrix.
   ///
   /// For real scalars a,b,c,d the Schur factorization
   ///
   ///
   ///        [ a  b ] = [ cs -sn ] [ A  B ] [ cs  sn ]
   ///        [ c  d ]   [ sn  cs ] [ C  D ] [-sn  cs ]
   /// is computed where either
   ///
   ///        1. C=0 so that A and D are real eigenvalues of the matrix, or
   ///        2. A=D and B*C<0, so that A +/- sqrt(B*C) are complex conjugate eigenvalues.
   /// @tparam real_t Real float point type.
   /// @param a Real scalar; matrix element a on input, matrix element A on exit.
   /// @param b Real scalar; matrix element b on input, matrix element B on exit.
   /// @param c Real scalar; matrix element c on input, matrix element C on exit.
   /// @param d Real scalar; matrix element d on input, matrix element D on exit.
   /// @param rt1r Real scalar; on exit the real part of first eigenvalue.
   /// @param rt1i Real scalar; on exit the imaginary part of first eigenvalue.
   /// @param rt2r Real scalar; on exit the real part of second eigenvalue.
   /// @param rt2i Real scalar; on exit the imaginary part of second eigenvalue.
   /// @param sn Real scalar; on exit the matrix element sn.
   /// @param cs Real scalar; on exit the matrix element cs.
   /// @ingroup AUX
   
   template<typename real_t>
   void LANV2(real_t &a,real_t &b,real_t &c,real_t &d,real_t &rt1r,real_t &rt1i,real_t &rt2r,real_t &rt2i,real_t &cs,real_t &sn)
   {
      using std::abs;
      using std::sqrt;
      using std::min;
      using std::max;
      using std::copysign;
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t half=0.5;
      const real_t multpl=4.0;
      const real_t eps=LAMCH<real_t>('P');
      if(c==zero)
      {
         cs=one;
         sn=zero;
      }
      else if(b==zero)
      {
         cs=zero;
         sn=one;
         real_t temp=d;
         d=a;
         a=temp;
         b=-c;
         c=zero;
      }
      else if(((a-d)==zero)&&(copysign(one,b)!=copysign(one,c)))
      {
         cs=one;
         sn=zero;
      }
      else
      {
         real_t temp=a-d;
         real_t p=half*temp;
         real_t bcmax=max(abs(b),abs(c));
         real_t bcmis=min(abs(b),abs(c))*copysign(one,b)*copysign(one,c);
         real_t scale=max(abs(p),bcmax);
         real_t z=(p/scale)*p+(bcmax/scale)*bcmis;
         if(z>=multpl*eps)
         {
            z=p+copysign(sqrt(scale)*sqrt(z),p);
            a=d+z;
            d=d-(bcmax/z)*bcmis;
            real_t tau=LAPY2(c,z);
            cs=z/tau;
            sn=c/tau;
            b=b-c;
            c=zero;
         }
         else
         {
            real_t sigma=b+c;
            real_t tau=LAPY2(sigma,temp);
            cs=sqrt(half*(one+abs(sigma)/tau));
            sn=-(p/(tau*cs))*copysign(one,sigma);
            real_t A=a*cs+b*sn;
            real_t B=-a*sn+b*cs;
            real_t C=c*cs+d*sn;
            real_t D=-c*sn+d*cs;
            a=A*cs+C*sn;
            b=B*cs+D*sn;
            c=-A*sn+C*cs;
            d=-B*sn+D*cs;
            temp=half*(a+d);
            a=temp;
            d=temp;
            if(c!=zero)
            {
               if(b!=zero)
               {
                  if(copysign(one,b)==copysign(one,c))
                  {
                     real_t sab=sqrt(abs(b));
                     real_t sac=sqrt(abs(c));
                     p=copysign(sab*sac,c);
                     tau=one/sqrt(abs(b+c));
                     a=temp+p;
                     d=temp-p;
                     b=b-c;
                     c=zero;
                     real_t cs1=sab*tau;
                     real_t sn1=sac*tau;
                     temp=cs*cs1-sn*sn1;
                     sn=cs*sn1+sn*cs1;
                     cs=temp;
                  }
               }
               else
               {
                  b=-c;
                  c=zero;
                  temp=cs;
                  cs=-sn;
                  sn=temp;
               }
            }
         }
      }
      rt1r=a;
      rt2r=d;
      if(c==zero)
      {
         rt1i=zero;
         rt2i=zero;
      }
      else
      {
         rt1i=sqrt(abs(b))*sqrt(abs(c));
         rt2i=-rt1i;
      }
   }
}

#endif
