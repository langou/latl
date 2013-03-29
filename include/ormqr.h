//
//  ormqr.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _ormqr_h
#define _ormqr_h

/// @file ormqr.h Applies orthogonal matrix Q to a matrix C using blocked algorithm.

#include <cctype>
#include <algorithm>
#include "larfb.h"
#include "latl.h"

namespace LATL
{
   /// @brief Applies orthogonal matrix Q to a matrix C.
   ///
   /// Overwrites the general real m-by-n matrix C with one of the following
   ///
   ///                     side='L'   side='R'
   ///        trans='N':   Q * C      C * Q
   ///        trans='T':   Q'* C      C * Q'
   /// where Q is a real orthogonal matrix defined as the product of k
   /// elementary reflectors
   ///
   ///        Q = H_1 H_2 . . . H_k
   /// as returned by LATL::GEQRF.  Q is of order m if side='L' and of order n if side='R'.
   /// @tparam real_t Floating point parameter.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @param side Specifies which side Q is to be applied.
   ///
   ///             'L': apply Q or Q' from the Left;
   ///             'R': apply Q or Q' from the Right.
   /// @param trans Specifies whether Q or Q' is applied.
   ///
   ///             'N':  No transpose, apply Q;
   ///             'T':  Transpose, apply Q'.
   /// @param m The number of rows of the matrix C.
   /// @param n The number of columns of the matrix C.
   /// @param k The number of elementary reflectors whose product defines the matrix Q.
   ///
   ///          If side='L', m >= k >= 0;
   ///          if side='R', n >= k >= 0.
   /// @param[in] A Real matrix containing the elementary reflectors H as returned by LATL::GEQRF.
   /// The i-th column must contain the vector which defines the elementary reflector H_i,
   /// for i = 0,2,...,k-1.  If side='L', A is m-by-k; if side='R', A is n-by-k.
   /// @param ldA The column length of the matrix A.  If side='L', ldA>=m; if side='R', ldA>=n.
   /// @param[in] tau Real vector of length k containing the scalar factors of the
   /// elementary reflectors as returned by LATL::GEQRF.
   /// @param[in,out] C Real m-by-n matrix.  On exit, C is overwritten by one of
   /// Q*C, Q'*C, C*Q*' or C*Q.
   /// @param ldC The column length the matrix C. ldC>=m.
   /// @param nb Block size.  0<nb<=k
   /// @ingroup COMP
   
   template<typename real_t>
   int ORMQR(char side, char trans, int_t m, int_t n int_t k, real_t *A, int_t ldA, real_t *tau, real_t *C, int_t ldC, int_t nb)
   {
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
      else if((k<0)||(k>((side=='L')?m:n)))
         return -5;
      else if(ldA<((side=='L')?m:n))
         return -7;
      else if(ldC<m)
         return -10;

      if((m==0)||(n==0)||(k==0))
         return 0;

      if((nb<1)||(nb>k))
         nb=k;

      real_t W=new real_t[n*nb+nb*nb];
      real_t *T=W+n*nb;
      const int_t ldT=nb;

      if((side=='L')&&(trans=='N'))
      {
         for(int_t i=0;i<k;i+=nb)
         {
            ib=min(nb,k-i);
            DLARFT('F','C',m-i,ib,A+i,ldA,tau,T,ldT);
            DLARFB('L','N','F','C',m-i,n-i,ib,A+i,ldA,T,ldT,C+i,ldC,W);
            A+=ib*ldA;
            tau+=ib;
         }
      }
      else if((side=='L')&&(trans=='T'))
      {

      }
      else if((side=='R')&&(trans=='N'))
      {

      }
      else // (side=='R')&&(trans=='T')
      {

      }
      delete [] W;
      return 0;
   }
}
#endif
