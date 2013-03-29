//
//  las2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _las2_h
#define _las2_h

/// @file las2.h Computes the singular values of a 2-by-2 triangular matrix.

#include <cmath>
#include <algorithm>

namespace LATL
{
   /// @brief Computes the singular values of a 2-by-2 triangular matrix.
   /// 
   /// Let f,g,h be real scalarsl the singular values of the matrix
   ///  
   ///       [  f   g  ]
   ///       [  0   h  ]
   /// are computed.
   /// @tparam real_t Floating point type.
   /// @param f Real scalar, the first element of the diagonal of the input matrix.
   /// @param g Real scalar, the superdiagaonal element of the input matrix.
   /// @param h Real scalar, the second element of the diagonal of the input matrix.
   /// @param ssmin On exit, the smaller singular value.
   /// @param ssmax On exit, the larger singular value.
   /// @ingroup AUX

   template<typename real_t>
   void LAS2(real_t f,real_t g,real_t h,real_t &ssmin,real_t &ssmax)
   {
      using std::abs;
      using std::max;
      using std::min;
      using std::sqrt;
      
      const real_t zero=0.0;
      const real_t one=1.0;
      const real_t two=2.0;
      
      real_t fa=abs(f);
      real_t ga=abs(g);
      real_t ha=abs(h);
      real_t fhmn=min(fa,ha);
      real_t fhmx=max(fa,ha);
      if(fhmn==zero)
      {
         ssmin=zero;
         if(fhmx==zero)
         {
            ssmax=ga;
         }
         else
         {
            real_t temp=min(fhmx,ga)/max(fhmx,ga);
            ssmax=max(fhmx,ga)*sqrt(one+temp*temp);
         }
      }
      else
      {
         if(ga<fhmx)
         {
            real_t as=one+fhmn/fhmx;
            real_t at=(fhmx-fhmn)/fhmx;
            real_t temp=ga/fhmx;
            real_t au=temp*temp;
            real_t c=two/(sqrt(as*as+au)+sqrt(at*at+au));
            ssmin=fhmn*c;
            ssmax=fhmx/c;
         }
         else
         {
            real_t au=fhmx/ga;
            if(au==zero)
            {
               ssmin=(fhmn*fhmx)/ga;
               ssmax=ga;
            }
            else
            {
               real_t as=one+fhmn/fhmx;
               real_t at=(fhmx-fhmn)/fhmx;
               real_t c=one/(sqrt(one+as*au*as*au)+sqrt(one+at*au*at*au));
               ssmin=(fhmn*c)*au;
               ssmin=ssmin+ssmin;
               ssmax=ga/(c+c);
            }
         }
      }
   }
}

#endif
