//
//  laneg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/20/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _laneg_h
#define _laneg_h

/// @file laneg.h Computes Sturm count of a matrix.

#include <algorithm>
#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Computes the Sturm count of a matrix.
   ///
   /// Computes the Sturm count, the number of negative pivots
   /// encountered while factoring tridiagonal 
   ///   
   ///        T - sigma I = L D L'.
   /// This implementation works directly on the factors without forming
   /// the tridiagonal matrix T.  The Sturm count is also the number of
   /// eigenvalues of T less than sigma.
   /// @return Sturm count.
   /// @tparam real_t Floating point type.
   /// @param n The order of the matrix.
   /// @param d Real vector of length n containing the diagonal matrix elements.
   /// @param lld Real vector of length n-1 containing the elements L(i)*L(i)*D(i).
   /// @param sigma Real scalar, specifies shift about in T - sigma I.
   /// @param r The twist index for the twisted factorization that is used for the negative count.
   /// @ingroup AUX

   template<typename real_t>
   int_t LANEG(int_t n, real_t *d, real_t *lld, real_t sigma, int_t r)
   {
      const real_t zero(0.0);
      const real_t one(1.0);
      const int_t blklen=128;
      int_t negcnt=0;

      // upper part: L D L' - sigma I = (L+)(D+)(L+')
      
      real_t t=-sigma;
      for(int_t bj=0;bj<r-1;bj+=blklen)
      {
         int_t neg1=0;
         real_t bsav=t;
         for(int_t j=bj;(j<bj+blklen)&&(j<r-1);j++)
         {
            real_t dplus=d[j]+t;
            if(dplus<zero)
               neg1++;
            real_t tmp=t/dplus;
            t=tmp*lld[j]-sigma;
         }
         if(std::isnan(t))
         {
            neg1=0;
            t=bsav;
            for(int_t j=bj;(j<bj+blklen)&&(j<r-1);j++)
            {
               real_t dplus=d[j]+t;
               if(dplus<zero)
                  neg1++;
               real_t tmp=t/dplus;
               if(std::isnan(tmp))
                  tmp=one;
               t=tmp*lld[j]-sigma;
            }
         }
         negcnt+=neg1;
      }
      
      // lower part:  L D L' - sigma I = (U-)(D-)(U-')
      
      real_t p=d[n-1]-sigma;
      for(int_t bj=n-2;bj>=r;bj-=blklen)
      {
         int_t neg2=0;
         real_t bsav=p;
         for(int_t j=bj;(j>=r)&&(j>bj-blklen);j--)
         {
            real_t dminus=lld[j]+p;
            if(dminus<zero)
               neg2++;
            real_t tmp=p/dminus;
            p=tmp*d[j]-sigma;
         }
         if(std::isnan(p))
         {
            neg2=0;
            p=bsav;
            for(int_t j=bj;(j>=r)&&(j>bj-blklen);j--)
            {
               real_t dminus=lld[j]+p;
               if(dminus<zero)
                  neg2++;
               real_t tmp=p/dminus;
               if(std::isnan(tmp))
                  tmp=one;
               p=tmp*d[j]-sigma;
            }
         }
         negcnt+=neg2;
      }
      
      // twist index t was shifted by sigma initially
      
      if((t+sigma)+p<zero)
         negcnt++;
      return negcnt;
   }
}

#endif
