//
//  ormr2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _ormr2_h
#define _ormr2_h

/// @file ormr2.h Applies orthogonal matrix Q to a matrix C.

#include <cctype>
#include <algorithm>
#include "larf.h"
#include "latl.h"

namespace LATL
{
   /// @brief Applies orthogonal matrix Q to a matrix C.
   /// @tparam real_t Floating point parameter.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @param side Specifies which side Q is to be applied.
   ///
   ///                 'L': apply Q or Q' from the Left;
   ///                 'R': apply Q or Q' from the Right.
   /// @param trans Specifies whether Q or Q' is applied.
   ///
   ///                 'N':  No transpose, apply Q;
   ///                 'T':  Transpose, apply Q'.
   /// @param m The number of rows of the matrix C.
   /// @param n The number of columns of the matrix C.
   /// @param k The number of elementary reflectors whose product defines the matrix Q.
   ///
   ///                 If side='L', m>=k>=0;
   ///                 if side='R', n>=k>=0.
   /// @param[in] A Real matrix containing the elementary reflectors H as
   /// returned by LATL::GEQR2 or LATL::GEQRF.
   ///
   ///                 If side='L', A is k-by-m;
   ///                 if side='R', A is k-by-n.
   /// @param ldA The column length of the matrix A.  ldA>=k.
   /// @param[in] tau Real vector of length k containing the scalar factors of the
   /// elementary reflectors as returned by LATL::GEQR2.
   /// @param[in,out] C Real m-by-n matrix.  On exit, C is replaced by one of the following:
   ///
   ///                 If side='L' & trans='N':  C <- Q * C
   ///                 If side='L' & trans='T':  C <- Q'* C
   ///                 If side='R' & trans='T':  C <- C * Q'
   ///                 If side='R' & trans='N':  C <- C * Q
   /// @param ldC The column length the matrix C. ldC>=m.
   /// @ingroup COMP
   
   template<typename real_t>
   int ORMR2(char side, char trans, int_t m, int_t n, int_t k, real_t *A, int_t ldA, real_t *tau, real_t *C, int_t ldC)
   {
      const real_t one(1.0);
      using std::toupper;
      using std::min;
      side=toupper(side);
      trans=toupper(trans);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((trans!='N')&&(trans!='T'))
         return -2;
      else if(m<0)
         return -3;
      else if(n<0)
         return -4;

      int_t q=(side=='L')?m:n;

      if((k<0)||(k>q))
         return -5;
      else if(ldA<q)
         return -7;
      else if(ldC<m)
         return -10;

      if((m==0)||(n==0)||(k==0))
         return 0;

      real_t *W=new real_t[q];
      real_t *v=new real_t[q];

      if((side=='L')&&(trans=='N'))
      {
         for(int_t i=0;i<k;i++)
         {
	   for(int_t j=0;j<=m-k+i;j++)
               v[j]=A[i+j*ldA];
            LARF('L',m-k+i+1,n,A+i,ldA,tau[i],C,ldC,W);
         }
      }
      else if((side=='L')&&(trans=='T'))
      {
         for(int_t i=k-1;i>=0;--i)
         {
            LARF('L',m-k+i+1,n,A+i,ldA,tau[i],C,ldC,W);
         }
      }
      else if((side=='R')&&(trans=='N'))
      {
         for(int_t i=0;i<k;i++)
         {
            LARF('R',m,n-k+i+1,A+i,ldA,tau[i],C,ldC,W);
         }
      }
      else // (side=='R')&&(trans=='T')
      {
         for(int_t i=k-1;i>=0;--i)
         {
            LARF('R',m,n-k+i+1,A+i,ldA,tau[i],C,ldC,W);
         }
      }
      delete [] v;
      delete [] W;
      return 0;
   }
}
#endif
