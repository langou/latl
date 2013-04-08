//
//  sytri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _sytri_h
#define _sytri_h

/// @file sytri.h Computes the inverse of a symmetric indefinite matrix.

#include "latl.h"
#include "dot.h"
#include "copy.h"
#include "swap.h"
#include "symv.h"
#include "syconv.h"
#include "trtri.h"
#include "trmm.h"
#include "syswapr.h"
#include <cmath>

namespace LATL
{
   /// @brief Computes the inverse of a symmetric indefinite matrix.
   ///
   /// A using the factorization A = U D U' or A = L D L' computed by
   /// LATL::SYTRF, the inverse of the symmetric indefinite matrix is returned in
   /// either the upper, U, or lower, L, part of the matrix A.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Real triangular matrix of order n.
   /// On entry, the block diagonal matrix D and the multipliers used to obtain
   /// the factor U or L as computed by LATL::SYTRF.  On exit, if upper 
   /// trianglar, A is overwritten with the upper triangle of the inverse of A;
   /// if lower trianglar, A is overwritten with the lower triangle of the
   /// inverse of A.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param ipiv is integer array with dimension n.  Details of the
   /// interchanges and the block structure of D as determined by LATL::SYTRF.
   /// @param bsdv Bool array size n.  On entry, contains the details of the block structure of D.
   /// @ingroup COMP
   
   template <typename real_t>
   int_t SYTRI(char uplo, int_t n, real_t *A, int_t ldA, int_t *ipiv, bool *bsdv)
   {
      using std::abs;
      using std::toupper;
      uplo=toupper(uplo);
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t t, Ak, Akp1, Akkp1, d, temp;

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;

      if(uplo=='U')
      {
         for(int_t i=n-1; i>=0; i--)
            if((bsdv[i] == 0) && (A[i+i*ldA]==zero))
               return 1;
      }
      else
      {
         for(int_t i=0; i<n; i++)
            if((bsdv[i] == 0) && (A[i+i*ldA]==zero))
               return 1;
      }

      real_t *work =new real_t[n];
      int_t kstep;

      if(uplo=='U')
      {
         int_t k = 1;
         while(k<=n)
         {
            if(bsdv[k-1] == 0)
            {
               A[(k-1)+(k-1)*ldA] = one / A[(k-1)+(k-1)*ldA];

               if(k>1)
               {
                  COPY<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( k-1, work, 1, A+(k-1)*ldA, 1);
               }
               kstep = 1;
            }
            else
            {
               t = abs( A[(k-1)+k*ldA] );
               Ak = A[(k-1)+(k-1)*ldA] / t;
               Akp1 = A[k+k*ldA] / t;
               Akkp1 = A[(k-1)+k*ldA] / t;
               d = t*( Ak*Akp1 - one );
               A[(k-1)+(k-1)*ldA] = Akp1 / d;
               A[k+k*ldA] = Ak / d;
               A[(k-1)+k*ldA] = -Akkp1 / d;

               if(k>1)
               {
                  COPY<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( k-1, work, 1, A+(k-1)*ldA, 1);
                  A[(k-1)+k*ldA] -= DOT( k-1, A+(k-1)*ldA, 1, A+k*ldA, 1);
                  COPY<real_t>( k-1, A+k*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+k*ldA, 1);
                  A[k+k*ldA] -= DOT( k-1, work, 1, A+k*ldA, 1);
               }
               kstep = 2;
            }

            int_t kp = ipiv[k-1];
            if(kp != k-1)
            {
               SWAP<real_t>( kp, A+(k-1)*ldA, 1, A+kp*ldA, 1);
               SWAP<real_t>( k-kp-2, A+(kp+1)+(k-1)*ldA, 1, A+kp+(kp+1)*ldA, ldA);
               temp = A[(k-1)+(k-1)*ldA];
               A[(k-1)+(k-1)*ldA] = A[kp+kp*ldA];
               A[kp+kp*ldA] = temp;
               if(kstep==2)
               {
                  temp = A[(k-1)+k*ldA];
                  A[(k-1)+k*ldA] = A[kp+k*ldA];
                  A[kp+k*ldA] = temp;
               }
            }
            k = k + kstep;
         }
      }
      else
      {
         int_t k = n;
         while(k>=1)
         {
            if(bsdv[k-1] == 0)
            {
               A[(k-1)+(k-1)*ldA] = one / A[(k-1)+(k-1)*ldA];

               if(k<n)
               {
                  COPY<real_t>( n-k, A+k+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( n-k, work, 1, A+k+(k-1)*ldA, 1);
               }
               kstep = 1;
            }
            else
            {
               t = abs( A[(k-1)+(k-2)*ldA] );
               Ak = A[(k-2)+(k-2)*ldA] / t;
               Akp1 = A[(k-1)+(k-1)*ldA] / t;
               Akkp1 = A[(k-1)+(k-2)*ldA] / t;
               d = t*( Ak*Akp1 - one );
               A[(k-2)+(k-2)*ldA] = Akp1 / d;
               A[(k-1)+(k-1)*ldA] = Ak / d;
               A[(k-1)+(k-2)*ldA] = -Akkp1 / d;

               if(k<n)
               {
                  COPY<real_t>( n-k, A+k+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( n-k, work, 1, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-2)*ldA] -= DOT( n-k, A+k+(k-1)*ldA, 1, A+k+(k-2)*ldA, 1);
                  COPY<real_t>( n-k, A+k+(k-2)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-2)*ldA, 1);
                  A[(k-2)+(k-2)*ldA] -= DOT( n-k, work, 1, A+k+(k-2)*ldA, 1);
               }
               kstep = 2;
            }

            int_t kp = ipiv[k-1];
            if(kp != k-1)
            {
               if(kp<n-1)
                  SWAP<real_t>( n-kp-1, A+kp+1+(k-1)*ldA,1,A+kp+1+kp*ldA,1);
               SWAP<real_t>( kp-k,A+k+(k-1)*ldA,1,A+kp+k*ldA,ldA);
               temp = A[(k-1)+(k-1)*ldA];
               A[(k-1)+(k-1)*ldA] = A[kp+kp*ldA];
               A[kp+kp*ldA] = temp;
               if(kstep==2)
               {
                  temp = A[(k-1)+(k-2)*ldA];
                  A[(k-1)+(k-2)*ldA] = A[kp+(k-2)*ldA];
                  A[kp+(k-2)*ldA] = temp;
               }
            }
            k = k - kstep;
         }
      }
      delete [] work;

      return 0;
   }

   /// @brief Computes the inverse of a symmetric indefinite matrix.
   ///
   /// A using the factorization A = U D U.' or A = L D L.' computed by
   /// LATL::SYTRF, the inverse of the symmetric indefinite matrix is returned in
   /// either the upper, U, or lower, L, part of the matrix A.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Complex triangular matrix of order n.
   /// On entry, the block diagonal matrix D and the multipliers used to obtain
   /// the factor U or L as computed by LATL::SYTRF.  On exit, if upper 
   /// trianglar, A is overwritten with the upper triangle of the inverse of A;
   /// if lower trianglar, A is overwritten with the lower triangle of the
   /// inverse of A.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param ipiv is integer array with dimension n.  Details of the
   /// interchanges and the block structure of D as determined by LATL::SYTRF.
   /// @param bsdv Bool array size n.  On entry, contains the details of the block structure of D.
   /// @ingroup COMP

   template <typename real_t>
   int_t SYTRI(char uplo, int_t n, complex<real_t> *A, int_t ldA, int_t *ipiv, bool *bsdv)
   {
      using std::toupper;
      uplo=toupper(uplo);
      const complex<real_t> zero(0.0);
      const complex<real_t> one(1.0);
      complex<real_t> t, Ak, Akp1, Akkp1, d, temp;

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;

      if(uplo=='U')
      {
         for(int_t i=n-1; i>=0; i--)
            if((bsdv[i] == 0) && (A[i+i*ldA]==zero))
               return 1;
      }
      else
      {
         for(int_t i=0; i<n; i++)
            if((bsdv[i] == 0) && (A[i+i*ldA]==zero))
               return 1;
      }


      complex<real_t> *work =new complex<real_t>[n];
      int_t kstep;

      if(uplo=='U')
      {
         int_t k = 1;
         while(k<=n)
         {
            if(bsdv[k-1] == 0)
            {
               A[(k-1)+(k-1)*ldA] = one / A[(k-1)+(k-1)*ldA];

               if(k>1)
               {
                  COPY<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( k-1, work, 1, A+(k-1)*ldA, 1);
               }
               kstep = 1;
            }
            else
            {
               t = A[(k-1)+k*ldA];
               Ak = A[(k-1)+(k-1)*ldA] / t;
               Akp1 = A[k+k*ldA] / t;
               Akkp1 = A[(k-1)+k*ldA] / t;
               d = t*( Ak*Akp1 - one );
               A[(k-1)+(k-1)*ldA] = Akp1 / d;
               A[k+k*ldA] = Ak / d;
               A[(k-1)+k*ldA] = -Akkp1 / d;

               if(k>1)
               {
                  COPY<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( k-1, work, 1, A+(k-1)*ldA, 1);
                  A[(k-1)+k*ldA] -= DOT( k-1, A+(k-1)*ldA, 1, A+k*ldA, 1);
                  COPY<real_t>( k-1, A+k*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+k*ldA, 1);
                  A[k+k*ldA] -= DOT( k-1, work, 1, A+k*ldA, 1);
               }
               kstep = 2;
            }

            int_t kp = ipiv[k-1];
            if(kp != k-1)
            {
               SWAP<real_t>( kp, A+(k-1)*ldA, 1, A+kp*ldA, 1);
               SWAP<real_t>( k-kp-2, A+(kp+1)+(k-1)*ldA, 1, A+kp+(kp+1)*ldA, ldA);
               temp = A[(k-1)+(k-1)*ldA];
               A[(k-1)+(k-1)*ldA] = A[kp+kp*ldA];
               A[kp+kp*ldA] = temp;
               if(kstep==2)
               {
                  temp = A[(k-1)+k*ldA];
                  A[(k-1)+k*ldA] = A[kp+k*ldA];
                  A[kp+k*ldA] = temp;
               }
            }
            k = k + kstep;
         }
      }
      else
      {
         int_t k = n;
         while(k>=1)
         {
            if(bsdv[k-1] == 0)
            {
               A[(k-1)+(k-1)*ldA] = one / A[(k-1)+(k-1)*ldA];

               if(k<n)
               {
                  COPY<real_t>( n-k, A+k+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( n-k, work, 1, A+k+(k-1)*ldA, 1);
               }
               kstep = 1;
            }
            else
            {
               t = A[(k-1)+(k-2)*ldA];
               Ak = A[(k-2)+(k-2)*ldA] / t;
               Akp1 = A[(k-1)+(k-1)*ldA] / t;
               Akkp1 = A[(k-1)+(k-2)*ldA] / t;
               d = t*( Ak*Akp1 - one );
               A[(k-2)+(k-2)*ldA] = Akp1 / d;
               A[(k-1)+(k-1)*ldA] = Ak / d;
               A[(k-1)+(k-2)*ldA] = -Akkp1 / d;

               if(k<n)
               {
                  COPY<real_t>( n-k, A+k+(k-1)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-1)*ldA] -= DOT( n-k, work, 1, A+k+(k-1)*ldA, 1);
                  A[(k-1)+(k-2)*ldA] -= DOT( n-k, A+k+(k-1)*ldA, 1, A+k+(k-2)*ldA, 1);
                  COPY<real_t>( n-k, A+k+(k-2)*ldA, 1, work, 1);
                  SYMV<real_t>( uplo, n-k, -one, A+k+k*ldA, ldA, work, 1, zero, A+k+(k-2)*ldA, 1);
                  A[(k-2)+(k-2)*ldA] -= DOT( n-k, work, 1, A+k+(k-2)*ldA, 1);
               }
               kstep = 2;
            }

            int_t kp = ipiv[k-1];
            if(kp != k-1)
            {
               if(kp<n-1)
                  SWAP<real_t>( n-kp-1, A+kp+1+(k-1)*ldA,1,A+kp+1+kp*ldA,1);
               SWAP<real_t>( kp-k,A+k+(k-1)*ldA,1,A+kp+k*ldA,ldA);
               temp = A[(k-1)+(k-1)*ldA];
               A[(k-1)+(k-1)*ldA] = A[kp+kp*ldA];
               A[kp+kp*ldA] = temp;
               if(kstep==2)
               {
                  temp = A[(k-1)+(k-2)*ldA];
                  A[(k-1)+(k-2)*ldA] = A[kp+(k-2)*ldA];
                  A[kp+(k-2)*ldA] = temp;
               }
            }
            k = k - kstep;
         }
      }

      delete [] work;
      
      return 0;
   }

}

#endif
