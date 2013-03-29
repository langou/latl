//
//  lasv2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lasv2_h
#define _lasv2_h

/// @file lasv2.h Computes the singular value decomposition of a 2-2 triangular matrix.

#include <cmath>
#include <limits>
#include "latl.h"

namespace LATL
{
   /// @brief Computes the singular value decomposition of a 2-by-2 triangular matrix.
   /// 
   /// Let f,g,h be real scalars, then the SVD of the matrix
   ///  
   ///       [  f   g  ]
   ///       [  0   h  ]
   /// is computed. On return, abs(ssmax) is the larger singular value, 
   /// abs(ssmin) is the smaller singular value, and (csl,snl) and (csr,snr) 
   /// are the left and right singular vectors for abs(ssmax), 
   /// giving the decomposition
   ///
   ///       [ csl  snl ] [  f   g  ] [ csr -snr ]  =  [ ssmax   0   ]
   ///       [-snl  csl ] [  0   h  ] [ snr  csr ]     [  0    ssmin ].
   /// @tparam real_t Floating point type.
   /// @param f Real scalar, the first element of the diagonal of the input matrix.
   /// @param g Real scalar, the superdiagaonal element of the input matrix.
   /// @param h Real scalar, the second element of the diagonal of the input matrix.
   /// @param ssmin On exit, the smaller in magnitude singular value.
   /// @param ssmax On exit, the larger in magnitude singular value.
   /// @param snr On exit, part of the right singular vector.
   /// @param csr On exit, part of the right singular vector.
   /// @param snl On exit, part of the left singular vector.
   /// @param csl On exit, part of the left singular vector.
   /// @ingroup AUX
   
   template<typename real_t>
   void LASV2(real_t f,real_t g,real_t h,real_t &ssmin,real_t &ssmax,real_t &snr,real_t &csr,real_t &snl,real_t &csl)
   {
      using std::sqrt;
      using std::abs;
      using std::copysign;
      using std::numeric_limits;
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t half=0.5;
      const real_t two=2.0;
      const real_t four=4.0;
      const real_t eps=numeric_limits<real_t>::epsilon();
      
      real_t ft=f;
      real_t fa=abs(ft);
      real_t ht=h;
      real_t ha=abs(h);
      real_t gt=g;
      real_t ga=abs(gt);      
      real_t clt,crt,slt,srt;
   
      int pmax=1;
      bool swap=(ha>fa);
      if(swap)
      {
         pmax=3;
         real_t temp=ft;
         ft=ht;
         ht=temp;
         temp=fa;
         fa=ha;
         ha=temp;
      }

      if(ga==zero)
      {
         ssmin=ha;
         ssmax=fa;
         clt=one;
         crt=one;
         slt=zero;
         srt=zero;
      }
      else
      {
         bool gasmal=1;
         if(ga>fa)
         {
            pmax=2;
            if((fa/ga)<eps)
            {
               gasmal=0;
               ssmax=ga;
               if(ga>one)
                  ssmin=fa/(ga/ha);
               else
                  ssmin=(fa/ga)*ha;
               clt=one;
               slt=ht/gt;
               srt=one;
               crt=ft/gt;
            }
         }
         if(gasmal)
         {
            real_t d=fa-ha;
            real_t l=(d==fa)?one:d/fa;
            real_t m=gt/ft;
            real_t t=two-l;
            real_t mm=m*m;
            real_t tt=t*t;
            real_t s=sqrt(tt+mm);
            real_t r=(l<zero)?abs(m):sqrt(l*l+mm);
            real_t a=half*(s+r);
            ssmin=ha/a;
            ssmax=fa*a;
            if(mm==zero)
            {
               if(l==zero)
                  t=copysign(two,ft)*copysign(one,gt);
               else
                  t=gt/copysign(d,ft)+m/t;
            }
            else
            {
               t=(m/(s+t)+m/(r+l))*(one+a);
            }
            l=sqrt(t*t+four);
            crt=two/l;
            srt=t/l;
            clt=(crt+srt*m)/a;
            slt=(ht/ft)*srt/a;
         }
      }
      if(swap)
      {
         csl=srt;
         snl=crt;
         csr=slt;
         snr=clt;
      }
      else
      {
         csl=clt;
         snl=slt;
         csr=crt;
         snr=srt;
      }
      real_t tsign=one;
      if(pmax==1)
         tsign=copysign(one,csr)*copysign(one,csl)*copysign(one,f);
      else if(pmax==2)
         tsign=copysign(one,snr)*copysign(one,csl)*copysign(one,g);
      else if(pmax==3)
         tsign=copysign(one,snr)*copysign(one,snl)*copysign(one,h);
         
      ssmax=copysign(ssmax,tsign);
      ssmin=copysign(ssmin,tsign*copysign(one,f)*copysign(one,h));
   }
}

#endif
