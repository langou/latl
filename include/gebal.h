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
   /// This involves, first, permuting A by a similarity transformation to isolate eigenvalues
   /// in the first 0 to ilo-1 and last ihi+1 to n-1 elements on the diagonal; and second,
   /// applying a diagonal similarity transformation to rows and columns ilo to ihi to make the
   /// rows and columns as close in norm as possible.  Both steps are optional.
   /// The permutations consist of row and column interchanges which put the matrix in the form
   ///
   ///                ( T1   X   Y  )
   ///       P' A P = (  0   B   Z  )
   ///                (  0   0   T2 )
   /// where T1 and T2 are upper triangular matrices whose eigenvalues lie
   /// along the diagonal.  The column indices ilo and ihi mark the starting
   /// and ending columns of the submatrix B. Balancing consists of applying
   /// a diagonal similarity transformation inv(D) * B * D to make the
   /// 1-norms of each row of B and its corresponding column nearly equal.
   /// The output matrix is
   ///
   ///                ( T1     X*D          Y    )
   ///                (  0  inv(D)*B*D  inv(D)*Z )
   ///                (  0      0           T2   ).
   /// Balancing may reduce the 1-norm of the matrix, and improve the
   /// accuracy of the computed eigenvalues and/or eigenvectors.
   /// See @cite Kressner05 pages 35-39, and @cite ParlettReinsch69.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith argument is invalid.
   /// @param job  Specifies the operations to be performed on A:
   ///
   ///        = 'N':  none; simply set ilo=0, ihi=n-1, D[i]=1 and P[i]=i for i=0..n-1;
   ///        = 'P':  permute only;
   ///        = 'S':  scale only;
   ///        = 'B':  both permute and scale.
   /// @param n The order of the matrix A.  n >= 0.
   /// @param[in,out] A Real n-by-n matrix.  On exit,  A is overwritten by the balanced matrix.
   /// @param ldA The leading dimension of the array A.  ldA >= n.
   /// @param[out] ilo Column index such that jth column of A is set to zero for j<ilo.
   /// @param[out] ihi Column index such that jth column of A is set to zero for j>ihi.
   /// @param[out] P Integer vector of length n containing the permutations applied to A.
   /// P[j] is the index of the row and column interchanged with row and column j.
   /// The order in which the interchanges are made is n-1 to ihi+1, then 0 to ilo-1.
   /// @param[out] D Real vector of length n containing the diagonal elements of the matrix D.
   /// D[j] is the scaling factor applied to row and column j.
   /// @ingroup COMP

   template<typename real_t>
   int GEBAL(char job, int_t n, real_t *A, int_t ldA, int_t &ilo, int_t &ihi, int_t *P, real_t *D)
   {
      using std::toupper;
      using std::abs;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      const real_t two(2.0);
      const real_t factor(0.95);
      const real_t beta(2.0);
      const real_t sfmin=beta*LAMCH<real_t>('S');
      const real_t sfmax=one/sfmin;

      using std::cout;
      using std::endl;
      
      job=toupper(job);
      if((job!='N')&&(job!='P')&&(job!='S')&&(job!='B'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;

      bool permute=0;
      bool scale=0;

      if((job=='P')||(job=='B'))
         permute=1;

      if((job=='S')||(job=='B'))
         scale=1;

      ilo=0;
      ihi=n-1;

      for(int_t i=0;i<n;i++)
      {
         D[i]=one;
         P[i]=i;
      }

      if(permute)
      {
         int_t i,j;
         bool swapped=1;
         while(swapped)
         {
            swapped=0;
            i=ilo;
            while((i<=ihi)&&(!swapped))
            {
               real_t s=zero;
               for(j=ilo;j<=ihi;j++)
               {
                  if(i==j) continue;
                  s+=abs(A[i+j*ldA]);
               }
               if(s==zero)
               {
                  SWAP<real_t>(n,A+i,ldA,A+ihi,ldA);
                  SWAP<real_t>(n,A+i*ldA,1,A+ihi*ldA,1);
                  P[ihi]=i;
                  ihi--;
                  swapped=1;
               }
               i++;
            }
            j=ilo;
            while((j<=ihi)&&(!swapped))
            {
               real_t s=zero;
               for(i=ilo;i<=ihi;i++)
               {
                  if(i==j) continue;
                  s+=abs(A[i+j*ldA]);
               }
               if(s==zero)
               {
                  SWAP<real_t>(n,A+ilo*ldA,1,A+j*ldA,1);
                  SWAP<real_t>(n,A+ilo,ldA,A+j,ldA);
                  P[ilo]=i;
                  ilo++;
                  swapped=1;
               }
               j++;
            }
         }
      }

      if(scale)
      {
         bool converged=0;
         while(!converged)
         {
            converged=1;
            for(int_t j=ilo;j<=ihi;j++)
            {
               real_t c=NRM2(ihi-ilo+1,A+ilo+j*ldA,1);
               real_t r=NRM2(ihi-ilo+1,A+j+ilo*ldA,ldA);
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
                        for(int_t i=ilo;i<=ihi;i++)
                           A[i+j*ldA]*=scal;
                        for(int_t i=ilo;i<=ihi;i++)
                           A[j+i*ldA]/=scal;
                     }
                  }
               }
            }
         }
      }
      return 0;
   }
}

#endif
