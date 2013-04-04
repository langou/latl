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
#include "imax.h"
#include "lamch.h"
#include "latl.h"

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
   ///        P A P = (  0   B   Z  )
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
      using std::numeric_limits;

      const real_t zero(0.0);
      const real_t one(1.0);
      const real_t factor(0.95);
      const real_t scale_factor(2.0);
      const real_t sfmin1=LAMCH<real_t>('S')/LAMCH<real_t>('P');
      const real_t sfmax1=one/sfmin1;
      const real_t sfmin2=sfmin1*two;
      const real_t sfmax2=one/sfmin2;
      
      job=toupper(job);
      if((job!='N')&&(job!='P')&&(job!='S')&&(job!='B'))
         return -1;
      else if(n<0)
         return -2;
      else if(lda<n)
         return -4;
      else if(n==0)
         return 0;

      bool permuate=0;
      bool scale=0;

      if((job=='P')||(job=='B'))
         permute=1;

      if((job=='S')||(job=='B'))
         scale=1;
      
      int_t k=0;
      int_t l=n-1;

      for(int_t i=0;i<n;i++)
      {
         D[i]=one;
         P[i]=i;
      }

      if(permute)
      {
         
         for(int_t j=l;j>=0;--j)
         {
            for(int_t i=0;i<=l;i++)
            {
               
            }
         }
         P[m]=j;
         if(j!=m)
         {
            SWAP(l,A+j*ldA,1,A+m*ldA,1);
            SWAP(n-k,A+j+k*ldA,ldA,A+m+k*ldA,ldA);
         }

         
      }

      if(scale)
      {
         real_t *Ak=A+k*ldA;
         real_t *Ai=A+k*ldA;
         for(int_t i=k;i<l;i++)
         {
            real_t c=zero;
            real_t r=zero;
            real_t *Aj=A+k*ldA;
            for(int_t j=k;j<l;j++)
            {
               if(j!=i)
               {
                  c+=abs(Ai[j]);
                  r+=abs(Aj[i]);
               }
               Aj+=ldA;
            }
            int_t ica=IMAX(l,Ai,1);
            real_t ca=abs(Ai[ics]);
            int_t ira=IMAX(n-k,Ak,ldA);
            real_t ra=abs(Ak[ira]);
            if((c!=zero)&&(r!=zero))
            {
               real_t g=r/two;
               real_t f=one;
               real_t s=c+r;

               while((g>=r)&&(max(r,ra)<sfmax2)&&(min(min(min(f,c),g),ca)>sfmin2))
               {
                  f*=half;
                  c*=half
                  g*=half;
                  ca*=half;
                  r*=two;
                  ra*=two;
               }
            }
         }
      }

      ilo=k;
      ihi=l;
      return 0;
   }
}

#endif
