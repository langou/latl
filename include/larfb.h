//
//  larfb.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 2/28/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _larfb_h
#define _larfb_h

/// @file larfb.h Applies a Householder block reflector to a matrix.

#include <cctype>
#include "gemm.h"
#include "trmm.h"
#include "lacpy.h"
#include "lapeq.h"
#include "latl.h"

namespace LATL
{

   /// @brief Applies a Householder block reflector to a real matrix.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the block relfector is applied to C from the left (side='L') or right (side='R').
   /// @param trans Specifies whether H (trans='N') or H' (trans='T') is applied to C.
   /// @param direct Specifies the storage direction of the Householder vectors in V: direct='F' for forward or direct='B' for backward.
   /// @param storeV Specifies whether the Householder vectors in V are stored in columns (storeV='C') or rows (storeV='R').
   /// @param m Number of rows of the matrix C.
   /// @param n Number of columns of the matrix C.
   /// @param k Order of the block reflector H.
   /// @param[in] V Matrix containing the Householder vectors.
   /// @param ldV Leading dimension of the matrix V.
   /// @param[in] T Triangular part of the block reflector of order k.
   /// @param ldT Leading dimension of matrix T.  ldT>=k
   /// @param[in,out] C Real m-by-n matrix.  On exit, C is overwritten with H*C, H'*C, C*H, or C*H'.
   /// @param ldC Leading dimension of matrix C.  ldC>=m
   /// @param W Workspace vector of length: k*n if side='L' or k*m if side='R'.
   /// @ingroup AUX
   
   template<typename real_t>
   void LARFB(char side, char trans, char direct, char storeV, int_t m, int_t n, int_t k, real_t *V, int_t ldV, real_t *T, int_t ldT, real_t *C, int_t ldC, real_t *W)
   {
      const real_t one(1.0);
      using std::toupper;
      side=toupper(side);
      trans=toupper(trans);
      direct=toupper(direct);
      storeV=toupper(storeV);

      if((m>0)&&(n>0))
      {
         if(storeV=='C')
         {
            if(direct=='F')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C,ldC,W,k);
                  TRMM('L','L','T','U',k,n,one,V,ldV,W,k);
                  if(m>k)
                     GEMM('T','N',k,n,m-k,one,V+k,ldV,C+k,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V+k,ldV,W,k,one,C+k,ldC);
                  TRMM('L','L','N','U',k,n,-one,V,ldV,W,k);
                  LAPEQ(k,n,C,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C,ldC,W,m);
                  TRMM('R','L','N','U',m,k,one,V,ldV,W,m);
                  if(n>k)
                     GEMM('N','N',m,k,n-k,one,C+k*ldC,ldC,V+k,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','T',m-k,n,k,-one,V+k,ldV,W,m,one,C+k,ldC);
                  TRMM('R','L','T','U',n,k,-one,V,ldV,W,m);
                  LAPEQ(n,k,C,ldC,one,W,m);
               }
            }
            else if(direct=='B')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C+m-k,ldC,W,k);
                  TRMM('L','L','T','U',k,n,one,V+m-k,ldV,W,k);
                  if(m>k)
                     GEMM('T','N',k,n,m-k,one,V,ldV,C,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V,ldV,W,k,one,C,ldC);
                  TRMM('L','L','N','U',k,n,-one,V+m-k,ldV,W,k);
                  LAPEQ(k,n,C+m-k,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C+(n-k)*ldC,ldC,W,m);
                  TRMM('R','L','N','U',m,k,one,V+n-k,ldV,W,m);
                  if(n>k)
                     GEMM('N','N',m,k,n-k,one,C,ldC,V,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','T',m-k,n,k,-one,V,ldV,W,m,one,C+k,ldC);
                  TRMM('R','L','T','U',n,k,-one,V+n-k,ldV,W,m);
                  LAPEQ(n,k,C+(n-k)*ldC,ldC,one,W,m);
               }
            }
         }
         else if(storeV=='R')
         {
            if(direct=='F')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C,ldC,W,k);
                  TRMM('L','U','N','U',k,n,one,V,ldV,W,k);
                  if(m>k)
                     GEMM('N','N',k,n,m-k,one,V+k*ldV,ldV,C+k,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('T','N',m-k,n,k,-one,V+k*ldV,ldV,W,k,one,C+k,ldC);
                  TRMM('L','U','T','U',k,n,-one,V,ldV,W,k);
                  LAPEQ(k,n,C,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C,ldC,W,m);
                  TRMM('R','U','T','U',m,k,one,V,ldV,W,m);
                  if(n>k)
                     GEMM('N','T',m,k,n-k,one,C+k*ldC,ldC,V+k*ldV,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','N',m,n-k,k,-one,W,m,V+k*ldV,ldV,one,C+k*ldC,ldC);
                  TRMM('R','U','N','U',m,k,-one,V,ldV,W,m);
                  LAPEQ(m,k,C,ldC,one,W,m);
               }
            }
            else if(direct=='B')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C+m-k,ldC,W,k);
                  TRMM('L','L','N','U',k,n,one,V+(m-k)*ldV,ldV,W,k);
                  if(m>k)
                     GEMM('N','N',k,n,m-k,one,V,ldV,C,ldC,one,W,k);
                  TRMM('L','L',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V,ldV,W,k,one,C,ldC);
                  TRMM('L','L','N','U',k,n,-one,V+(m-k)*ldV,ldV,W,k);
                  LAPEQ(k,n,C+m-k,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C+(n-k)*ldC,ldC,W,m);
                  TRMM('R','L','T','U',m,k,one,V+(n-k)*ldV,ldV,W,m);
                  if(n>k)
                     GEMM('N','T',m,k,n-k,one,C,ldC,V,ldV,one,W,m);
                  TRMM('R','L',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','N',m,n-k,k,-one,W,m,V,ldV,one,C,ldC);
                  TRMM('R','L','N','U',m,k,-one,V+(n-k)*ldV,ldV,W,m);
                  LAPEQ(n,k,C+(n-k)*ldC,ldC,one,W,m);
               }
            }
         }
      }
   }
   
   /// @brief Applies a Householder block reflector to a complex matrix.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the block relfector is applied to C from the left (side='L') or right (side='R').
   /// @param trans Specifies whether H (trans='N') or H' (trans='C') is applied to C.
   /// @param direct Specifies the storage direction of the Householder vectors in V: direct='F' for forward or direct='B' for backward.
   /// @param storeV Specifies whether the Householder vectors in V are stored in columns (storeV='C') or rows (storeV='R').
   /// @param m Number of rows of the matrix C.
   /// @param n Number of columns of the matrix C.
   /// @param k Order of the block reflector H.
   /// @param[in] V Matrix containing the Householder vectors.
   /// @param ldV Leading dimension of the matrix V.
   /// @param[in] T Triangular part of the block reflector of order k.
   /// @param ldT Leading dimension of matrix T.  ldT>=k
   /// @param[in,out] C Complex m-by-n matrix.  On exit, C is overwritten with H*C, H'*C, C*H, or C*H'.
   /// @param ldC Leading dimension of matrix C.  ldC>=m
   /// @param W Workspace vector of length: k*n if side='L' or k*m if side='R'.
   /// @ingroup AUX

   template<typename real_t>
   void LARFB(char side, char trans, char direct, char storeV, int_t m, int_t n, int_t k, complex<real_t> *V, int_t ldV, complex<real_t> *T, int_t ldT, complex<real_t> *C, int_t ldC, complex<real_t> *W)
   {
      const complex<real_t> one(1.0);
      using std::toupper;
      side=toupper(side);
      trans=toupper(trans);
      direct=toupper(direct);
      storeV=toupper(storeV);

      if((m>0)&&(n>0))
      {
         if(storeV=='C')
         {
            if(direct=='F')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C,ldC,W,k);
                  TRMM('L','L','C','U',k,n,one,V,ldV,W,k);
                  if(m>k)
                     GEMM('C','N',k,n,m-k,one,V+k,ldV,C+k,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V+k,ldV,W,k,one,C+k,ldC);
                  TRMM('L','L','N','U',k,n,-one,V,ldV,W,k);
                  LAPEQ(k,n,C,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C,ldC,W,m);
                  TRMM('R','L','N','U',m,k,one,V,ldV,W,m);
                  if(n>k)
                     GEMM('N','N',m,k,n-k,one,C+k*ldC,ldC,V+k,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','C',m-k,n,k,-one,V+k,ldV,W,m,one,C+k,ldC);
                  TRMM('R','L','C','U',n,k,-one,V,ldV,W,m);
                  LAPEQ(n,k,C,ldC,one,W,m);
               }
            }
            else if(direct=='B')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C+m-k,ldC,W,k);
                  TRMM('L','L','C','U',k,n,one,V+m-k,ldV,W,k);
                  if(m>k)
                     GEMM('C','N',k,n,m-k,one,V,ldV,C,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V,ldV,W,k,one,C,ldC);
                  TRMM('L','L','N','U',k,n,-one,V+m-k,ldV,W,k);
                  LAPEQ(k,n,C+m-k,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C+(n-k)*ldC,ldC,W,m);
                  TRMM('R','L','N','U',m,k,one,V+n-k,ldV,W,m);
                  if(n>k)
                     GEMM('N','N',m,k,n-k,one,C,ldC,V,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','C',m-k,n,k,-one,V,ldV,W,m,one,C+k,ldC);
                  TRMM('R','L','C','U',n,k,-one,V+n-k,ldV,W,m);
                  LAPEQ(n,k,C+(n-k)*ldC,ldC,one,W,m);
               }
            }
         }
         else if(storeV=='R')
         {
            if(direct=='F')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C,ldC,W,k);
                  TRMM('L','U','N','U',k,n,one,V,ldV,W,k);
                  if(m>k)
                     GEMM('N','N',k,n,m-k,one,V+k*ldV,ldV,C+k,ldC,one,W,k);
                  TRMM('L','U',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('C','N',m-k,n,k,-one,V+k*ldV,ldV,W,k,one,C+k,ldC);
                  TRMM('L','U','C','U',k,n,-one,V,ldV,W,k);
                  LAPEQ(k,n,C,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C,ldC,W,m);
                  TRMM('R','U','C','U',m,k,one,V,ldV,W,m);
                  if(n>k)
                     GEMM('N','C',m,k,n-k,one,C+k*ldC,ldC,V+k*ldV,ldV,one,W,m);
                  TRMM('R','U',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','N',m,n-k,k,-one,W,m,V+k*ldV,ldV,one,C+k*ldC,ldC);
                  TRMM('R','U','N','U',m,k,-one,V,ldV,W,m);
                  LAPEQ(m,k,C,ldC,one,W,m);
               }
            }
            else if(direct=='B')
            {
               if(side=='L')
               {
                  LACPY('A',k,n,C+m-k,ldC,W,k);
                  TRMM('L','L','N','U',k,n,one,V+(m-k)*ldV,ldV,W,k);
                  if(m>k)
                     GEMM('N','N',k,n,m-k,one,V,ldV,C,ldC,one,W,k);
                  TRMM('L','L',trans,'N',k,n,one,T,ldT,W,k);
                  if(m>k)
                     GEMM('N','N',m-k,n,k,-one,V,ldV,W,k,one,C,ldC);
                  TRMM('L','L','N','U',k,n,-one,V+(m-k)*ldV,ldV,W,k);
                  LAPEQ(k,n,C+m-k,ldC,one,W,k);
               }
               else if(side=='R')
               {
                  LACPY('A',m,k,C+(n-k)*ldC,ldC,W,m);
                  TRMM('R','L','C','U',m,k,one,V+(n-k)*ldV,ldV,W,m);
                  if(n>k)
                     GEMM('N','C',m,k,n-k,one,C,ldC,V,ldV,one,W,m);
                  TRMM('R','L',trans,'N',m,k,one,T,ldT,W,m);
                  if(n>k)
                     GEMM('N','N',m,n-k,k,-one,W,m,V,ldV,one,C,ldC);
                  TRMM('R','L','N','U',m,k,-one,V+(n-k)*ldV,ldV,W,m);
                  LAPEQ(n,k,C+(n-k)*ldC,ldC,one,W,m);
               }
            }
         }
      }
   }
   
}

#endif
