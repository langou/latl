//
//  unmr2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _unmr2_h
#define _unmr2_h

/// @file unmr2.h Applies unitary matrix Q to a matrix C.

#include <cctype>
#include <algorithm>
#include "larf.h"
#include "latl.h"

namespace LATL
{
   /// @brief Applies unitary matrix Q to a matrix C.
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
   ///                 'C':  Conjugate transpose, apply Q'.
   /// @param m The number of rows of the matrix C.
   /// @param n The number of columns of the matrix C.
   /// @param k The number of elementary reflectors whose product defines the matrix Q.
   ///
   ///                 If side='L', m>=k>=0;
   ///                 if side='R', n>=k>=0.
   /// @param[in] A Complex matrix containing the elementary reflectors H as
   /// returned by LATL::GEQR2 or LATL::GEQRF.
   ///
   ///                 If side='L', A is k-by-m;
   ///                 if side='R', A is k-by-n.
   /// @param ldA The column length of the matrix A.  ldA>=k.
   /// @param[in] tau Complex vector of length k containing the scalar factors of the
   /// elementary reflectors as returned by LATL::GEQR2.
   /// @param[in,out] C Complex m-by-n matrix.  On exit, C is replaced by one of the following:
   ///
   ///                 If side='L' & trans='N':  C <- Q * C
   ///                 If side='L' & trans='C':  C <- Q'* C
   ///                 If side='R' & trans='C':  C <- C * Q'
   ///                 If side='R' & trans='N':  C <- C * Q
   /// @param ldC The column length the matrix C. ldC>=m.
   /// @param W Workspace vector (optional).
   ///
   ///                 If side='L', length of W is m;
   ///                 if side='R', length of W is n.
   /// If not provided, workspace is managed internally.
   /// @ingroup COMP
   
   template<typename real_t>
   int unmr2(char side, char trans, int_t m, int_t n, int_t k, complex<real_t> *A, int_t ldA, complex<real_t> *tau, complex<real_t> *C, int_t ldC, complex<real_t> *W=NULL)
   {
      const complex<real_t> one(1.0);
      using std::toupper;
      using std::min;
      using std::conj;
      side=toupper(side);
      trans=toupper(trans);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((trans!='N')&&(trans!='C'))
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

      bool allocate=(W==NULL)?1:0;
      if(allocate)
         W=new complex<real_t>[(side=='L')?m:n];

      if((side=='L')&&(trans=='N'))
      {
         complex<real_t> *B=A+(m-k)*ldA;
         for(int_t i=0;i<k;i++)
         {
            complex<real_t> alpha=B[i];
            B[i]=one;
            LARF('L',m-k+i+1,n,A+i,ldA,tau[i],C,ldC,W);
            B[i]=alpha;
            B+=ldA;
         }
      }
      else if((side=='L')&&(trans=='C'))
      {
         complex<real_t> *B=A+m*ldA;
         for(int_t i=k-1;i>=0;--i)
         {
            B-=ldA;
            complex<real_t> alpha=B[i];
            B[i]=one;
            LARF('L',m-k+i+1,n,A+i,ldA,conj(tau[i]),C,ldC,W);
            B[i]=alpha;
         }
      }
      else if((side=='R')&&(trans=='N'))
      {
         complex<real_t> *B=A+(n-k)*ldA;
         for(int_t i=0;i<k;i++)
         {
            complex<real_t> alpha=B[i];
            B[i]=one;
            LARF('R',m,n-k+i+1,A+i,ldA,tau[i],C,ldC,W);
            B[i]=alpha;
            B+=ldA;
         }
      }
      else // (side=='R')&&(trans=='C')
      {
         complex<real_t> *B=A+n*ldA;
         for(int_t i=k-1;i>=0;--i)
         {
            B-=ldA;
            LACGV(m,n-k+i+1,A+i,ldA);
            complex<real_t> alpha=B[i];
            B[i]=one;
            LARF('R',m,n-k+i+1,A+i,ldA,conj(tau[i]),C,ldC,W);
            B[i]=alpha;
         }
      }
      if(allocate)
         delete [] W;
      return 0;
   }
}
#endif
