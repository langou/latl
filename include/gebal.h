//
//  gebal.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/5/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _gebal_h
#define _gebal_h

/// @file gebal.h Balances a general matrix.

#include <cctype>
#include <cmath>
#include "swap.h"
#include "lamch.h"
#include "nrm2.h"
#include "latl.h"

#include <iostream>

namespace LATL
{
   /// @brief Balances a general real matrix A.
   ///
   /// Balancing consists of applying a diagonal similarity transformation
   ///
   ///        B = D^(-1) * A * D
   /// to make the norms of each row of B and its corresponding column nearly equal.
   /// Balancing may reduce the norm of the matrix, and improve the
   /// accuracy of the computed eigenvalues and/or eigenvectors.
   /// See @cite Kressner05 pages 35-39, and @cite ParlettReinsch69.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith argument is invalid.
   /// @param n The order of the matrix A.  n >= 0.
   /// @param[in,out] A Real n-by-n matrix.  On exit,  A is overwritten by the balanced matrix B.
   /// @param ldA The leading dimension of the array A.  ldA >= n.
   /// @param[out] D Real vector of length n containing the diagonal elements of the matrix D.
   /// @ingroup COMP

   template<typename real_t>
   int GEBAL(int_t n, real_t *A, int_t ldA, real_t *D)
   {
      using std::abs;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      const real_t two(2.0);
      const real_t factor(0.95);
      const real_t beta(2.0);
      const real_t sfmin=beta*LAMCH<real_t>('S');
      const real_t sfmax=one/sfmin;

      if(n<0)
         return -1;
      else if(ldA<n)
         return -3;
      else if(n==0)
         return 0;

      for(int_t i=0;i<n;i++)
         D[i]=one;

      bool converged=0;
      while(!converged)
      {
         converged=1;
         for(int_t j=0;j<n;j++)
         {
            real_t c=NRM2(n,A+j*ldA,1);
            real_t r=NRM2(n,A+j,ldA);
            real_t scal=one;
            if(isnan(c+r))
               return -3;
            if((c>zero)&&(r>zero))
            {
               real_t s=c*c+r*r;
               while((c<r/beta)&&(c<sfmax)&&(scal<sfmax)&&(r>sfmin))
               {
                  c*=beta;
                  r/=beta;
                  scal*=beta;
               }
               while((c>r*beta)&&(c>sfmin)&&(scal>sfmin)&&(r<sfmax))
               {
                  c/=beta;
                  r*=beta;
                  scal/=beta;
               }
               if((c*c+r*r)<factor*s)
               {
                  if(!((scal<one)&&(D[j]<one)&&(D[j]*scal<sfmin))&&!((scal>one)&&(D[j]>one)&&(D[j]>sfmax/scal)))
                  {
                     converged=0;
                     D[j]*=scal;
                     for(int_t i=0;i<n;i++)
                        A[i+j*ldA]*=scal;
                     for(int_t i=0;i<n;i++)
                        A[j+i*ldA]/=scal;
                  }
               }
            }
         }
      }
      return 0;
   }
}

#endif
