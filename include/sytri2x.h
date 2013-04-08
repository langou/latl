//
//  sytri2x.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 4/03/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _sytri2x_h
#define _sytri2x_h

/// @file sytri2x.h Computes the inverse of a symmetric indefinite matrix.

#include "latl.h"
#include "dot.h"
#include "copy.h"
#include "swap.h"
#include "symv.h"
#include "syconv.h"
#include "trtri.h"
#include "sytri.h"
#include "trmm.h"
#include "syswapr.h"
#include <cmath>

namespace LATL
{
   /// @brief Computes the inverse of a real symmetric indefinite matrix via
   /// blocked routine.
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
   /// @param bsdv is a boolean array with dimension n.  Details of the 
   /// interchanges and the block structure of D as determined by LATL::SYTRF.
   /// @param nb Block size.
   /// @ingroup COMP

   template <typename real_t>
   int_t SYTRI2X(char uplo, int_t n, real_t *A, int_t ldA, int_t *ipiv, bool *bsdv, int_t nb)
   {
      using std::abs;
      using std::toupper;
      int_t i, u11, ip, 
            k, cut, 
            nnb, count, j, 
            invd, ldwork;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t Ak, Akkp1, Akp1, d, t,
             u01_i_j, u01_ip1_j,
             u11_i_j, u11_ip1_j;

      uplo = toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;
      else if((nb<=1) || (nb>n))
      {
         LATL::SYTRI<real_t>( uplo, n, A, ldA, ipiv, bsdv );
         return 0;
      }

      ldwork = n + nb + 1;
      real_t *work = new real_t[ldwork*(nb+3)];
      SYCONV<real_t>( uplo, 'C', n, A, ldA, ipiv, bsdv, work );

      real_t *Aii;
      if(uplo == 'U')
      {
         Aii = A+(n-1)*(1+ldA);
         for(i=n-1;i>=0;i--)
         {
            if((bsdv[i] == 0) && (Aii[0] == zero))
               return -1;
            Aii--;
         }
      }
      else
      {
         Aii = A;
         for(i=0;i<n;i++)
         {
            if((bsdv[i] == 0) && (Aii[0] == zero))
               return -1;
            Aii++;
         }
      }

      u11 = n-1;
      invd = nb + 2 - 1;

      if(uplo == 'U')
      {
         TRTRI<real_t>( uplo, 'U', n, A, ldA, nb );

         k = 0;
         while(k<n)
         {
            if(bsdv[k] == 0)
            {
               work[k+invd*ldwork] = one / A[k+k*ldA];
               work[k+(invd+1)*ldwork] = zero;
               k++;
            }
            else
            {
               t = work[k+1];
               Ak = A[k+k*ldA] / t;
               Akp1 = A[k+1+(k+1)*ldA] / t;
               Akkp1 = work[k+1] / t;
               d = t * ( Ak * Akp1 - one );
               work[k+invd*ldwork] = Akp1 / d;
               work[k+1+(invd+1)*ldwork] = Ak / d;
               work[k+(invd+1)*ldwork] = -Akkp1 / d;
               work[k+1+invd*ldwork] = -Akkp1 / d;
               k+=2;
            }
         }

         cut = n;
         while(cut > 0)
         {
            nnb = nb;
            if( cut <= nnb)
               nnb = cut;
            else
            {
               count = 0;
               for(i=cut-nnb;i<cut;i++)
                  if (bsdv[i] == 1)
                     count++;
               if(count%2 == 1)
                  nnb++;
            }

            cut-=nnb;

            for(i=0;i<cut;i++)
               for(j=0;j<nnb;j++)
                  work[i+j*ldwork] = A[i+(cut+j)*ldA];

            for(i=0;i<nnb;i++)
            {
               work[u11+i+1+i*ldwork] = one;
               for(j=0;j<i;j++)
                  work[u11+i+1+j*ldwork] = zero;
               for(j=i+1;j<nnb;j++)
                  work[u11+i+1+j*ldwork] = A[cut+i+(cut+j)*ldA];
            }

            i = 0;
            while(i<cut)
            {
               if(bsdv[i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[i+j*ldwork] = work[i+invd*ldwork]*work[i+j*ldwork];
                  i++;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u01_i_j = work[i+j*ldwork];
                     u01_ip1_j = work[i+1+j*ldwork];
                     work[i+j*ldwork] = work[i+invd*ldwork] * u01_i_j + work[i+(invd+1)*ldwork] * u01_ip1_j;
                     work[i+1+j*ldwork] = work[i+1+invd*ldwork] * u01_i_j + work[i+1+(invd+1)*ldwork] * u01_ip1_j;
                  }
                  i+=2;
               }
            }

            i = 0;
            while(i<nnb)
            {
               if(bsdv[cut+i] == 0)
               {
                  for(j=i;j<nnb;j++)
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork];
                  i++;
               }
               else
               {
                  for(j=i;j<nnb;j++)
                  {
                     u11_i_j = work[u11+i+1+j*ldwork];
                     u11_ip1_j = work[u11+i+2+j*ldwork];
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork] + work[cut+i+(invd+1)*ldwork] * work[u11+i+2+j*ldwork];
                     work[u11+i+2+j*ldwork] = work[cut+i+1+invd*ldwork] * u11_i_j + work[cut+i+1+(invd+1)*ldwork] * u11_ip1_j;
                  }
                  i+=2;
               }
            }

            TRMM<real_t>( 'L', 'U', 'T', 'U', nnb, nnb, one, A+cut+cut*ldA, ldA, work+u11+1, n+nb+1 );
            for(i=0;i<nnb;i++)
               for(j=i;j<nnb;j++)
                  A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];

            GEMM<real_t>( 'T', 'N', nnb, nnb, cut, one, A+cut*ldA, ldA, work, n+nb+1, zero, work+u11+1, n+nb+1 );
            for(i=0;i<nnb;i++)
               for(j=i;j<nnb;j++)
                  A[cut+i+(cut+j)*ldA] = A[cut+i+(cut+j)*ldA] + work[u11+i+1+j*ldwork];

            TRMM<real_t>( 'L', 'U', 'T', 'U', cut, nnb, one, A, ldA, work, n+nb+1 );

            for(i=0;i<cut;i++)
               for(j=0;j<nnb;j++)
                  A[i+(cut+j)*ldA] = work[i+j*ldwork];
         }

         i = 0;
         while(i<n)
         {
            ip = ipiv[i];
            if(bsdv[i] == 0)
            {
               if(i<ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i>ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
            }
            else
            {
               i++;
               if((i-1) < ip)
                  SYSWAPR( uplo, n, A, ldA, i-1, ip );
               if((i-1) > ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i-1 );
            }
            i++;
         }
      }
      else
      {

         TRTRI<real_t>( 'L', 'U', n, A, ldA, nb );

         k = n-1;
         while(k>=0)
         {
            if(bsdv[k] == 0)
            {
               work[k+invd*ldwork] = one / A[k+k*ldA];
               work[k+(invd+1)*ldwork] = zero;
               k--;
            }
            else
            {
               t = work[k-1];
               Ak = A[k-1+(k-1)*ldA] / t;
               Akp1 = A[k+k*ldA] / t;
               Akkp1 = work[k-1] / t;
               d = t * ( Ak * Akp1 - one );
               work[k-1+invd*ldwork] = Akp1 / d;
               work[k+invd*ldwork] = Ak / d;
               work[k+(invd+1)*ldwork] = -Akkp1 / d;
               work[k-1+(invd+1)*ldwork] = -Akkp1 / d;
               k-=2;
            }
         }

         cut = 0;
         while(cut < n)
         {
            nnb = nb;
            if(cut+nnb > n)
               nnb = n - cut;
            else
            {
               count = 0;
               for(i=cut;i<cut+nnb;i++)
                  if (bsdv[i] == 1)
                     count++;
               if(count%2 == 1)
                  nnb++;
            }

            for(i=0;i<n-cut-nnb;i++)
               for(j=0;j<nnb;j++)
                  work[i+j*ldwork] = A[cut+nnb+i+(cut+j)*ldA];

            for(i=0;i<nnb;i++)
            {
               work[u11+i+1+i*ldwork] = one;
               for(j=i+1;j<nnb;j++)
                  work[u11+i+1+j*ldwork] = zero;
               for(j=0;j<=i-1;j++)
                  work[u11+i+1+j*ldwork] = A[cut+i+(cut+j)*ldA];
            }

            i = n - 1 - cut - nnb;
            while(i>=0)
            {
               if(bsdv[cut+nnb+i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[i+j*ldwork] = work[cut+nnb+i+invd*ldwork]*work[i+j*ldwork];
                  i--;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u01_i_j = work[i+j*ldwork];
                     u01_ip1_j = work[i-1+j*ldwork];
                     work[i+j*ldwork] = work[cut+nnb+i+invd*ldwork] * u01_i_j + work[cut+nnb+i+(invd+1)*ldwork] * u01_ip1_j;
                     work[i-1+j*ldwork] = work[cut+nnb+i-1+(invd+1)*ldwork] * u01_i_j + work[cut+nnb+i-1+invd*ldwork] * u01_ip1_j;
                  }
                  i-=2;
               }
            }

            i = nnb - 1;
            while(i>=0)
            {
               if(bsdv[cut+i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork];
                  i--;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u11_i_j = work[u11+i+1+j*ldwork];
                     u11_ip1_j = work[u11+i+j*ldwork];
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork] + work[cut+i+(invd+1)*ldwork] * u11_ip1_j;
                     work[u11+i+j*ldwork] = work[cut+i-1+(invd+1)*ldwork] * u11_i_j + work[cut+i-1+invd*ldwork] * u11_ip1_j;
                  }
                  i-=2;
               }
            }

            TRMM<real_t>( 'L', 'L', 'T', 'U', nnb, nnb, one, A+cut+cut*ldA, ldA, work+u11+1, n+nb+1 );

            for(i=0;i<nnb;i++)
               for(j=0;j<=i;j++)
                  A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];

            if((cut+nnb) < n)
            {
               GEMM<real_t>( 'T', 'N', nnb, nnb, n-nnb-cut, one, A+cut+nnb+cut*ldA, ldA, work, n+nb+1, zero, work+u11+1, n+nb+1 );

               for(i=0;i<nnb;i++)
                  for(j=0;j<=i;j++)
                     A[cut+i+(cut+j)*ldA] = A[cut+i+(cut+j)*ldA] + work[u11+i+1+j*ldwork];

               TRMM<real_t>( 'L', 'L', 'T', 'U', n-nnb-cut, nnb, one, A+cut+nnb+(cut+nnb)*ldA, ldA, work, n+nb+1 );

               for(i=0;i<n-cut-nnb;i++)
                  for(j=0;j<nnb;j++)
                     A[cut+nnb+i+(cut+j)*ldA] = work[i+j*ldwork];
            }
            else
            {
               for(i=0;i<nnb;i++)
                  for(j=0;j<i;j++)
                     A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];
            }

            cut += nnb;
         }

         i = n - 1;
         while(i>=0)
         {
            ip = ipiv[i];
            if(bsdv[i] == 0)
            {
               if(i<ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i>ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
            }
            else
            {
               if(i < ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i > ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
               i--;
            }
            i--;
         }
      }

      delete [] work;
      return 0;
   }

   /// @brief Computes the inverse of a complex symmetric indefinite matrix.
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
   /// @param A Complex triangular matrix of order n.
   /// On entry, the block diagonal matrix D and the multipliers used to obtain
   /// the factor U or L as computed by LATL::SYTRF.  On exit, if upper 
   /// trianglar, A is overwritten with the upper triangle of the inverse of A;
   /// if lower trianglar, A is overwritten with the lower triangle of the
   /// inverse of A.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param ipiv is integer array with dimension n.  Details of the
   /// interchanges and the block structure of D as determined by LATL::SYTRF.
   /// @param bsdv is a boolean array with dimension n.  Details of the 
   /// interchanges and the block structure of D as determined by LATL::SYTRF.
   /// @param nb Block size.
   /// @ingroup COMP

   template <typename real_t>
   int_t SYTRI2X(char uplo, int_t n, complex<real_t> *A, int_t ldA, int_t *ipiv, bool *bsdv, int_t nb)
   {
      using std::abs;
      using std::toupper;
      int_t i, u11, ip, 
            k, cut, 
            nnb, count, j, 
            invd, ldwork;
      const complex<real_t> zero(0.0);
      const complex<real_t> one(1.0);
      complex<real_t> Ak, Akkp1, Akp1, d, t,
             u01_i_j, u01_ip1_j,
             u11_i_j, u11_ip1_j;

      uplo = toupper(uplo);
      ldwork = n + nb + 1;

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;
      else if((nb<=1) || (nb>n))
      {
         LATL::SYTRI<real_t>( uplo, n, A, ldA, ipiv, bsdv );
         return 0;
      }

      complex<real_t> *work = new complex<real_t>[ldwork*(nb+3)];
      SYCONV<real_t>( uplo, 'C', n, A, ldA, ipiv, bsdv, work );

      complex<real_t> *Aii;
      if(uplo == 'U')
      {
         Aii = A+(n-1)*(1+ldA);
         for(i=n-1;i>=0;i--)
         {
            if((bsdv[i] == 0) && (Aii[0] == zero))
               return -1;
            Aii--;
         }
      }
      else
      {
         Aii = A;
         for(i=0;i<n;i++)
         {
            if((bsdv[i] == 0) && (Aii[0] == zero))
               return -1;
            Aii++;
         }
      }

      u11 = n-1;
      invd = nb + 2 - 1;

      if(uplo == 'U')
      {
         TRTRI<real_t>( uplo, 'U', n, A, ldA, nb );

         k = 0;
         while(k<n)
         {
            if(bsdv[k] == 0)
            {
               work[k+invd*ldwork] = one / A[k+k*ldA];
               work[k+(invd+1)*ldwork] = zero;
               k++;
            }
            else
            {
               t = work[k+1];
               Ak = A[k+k*ldA] / t;
               Akp1 = A[k+1+(k+1)*ldA] / t;
               Akkp1 = work[k+1] / t;
               d = t * ( Ak * Akp1 - one );
               work[k+invd*ldwork] = Akp1 / d;
               work[k+1+(invd+1)*ldwork] = Ak / d;
               work[k+(invd+1)*ldwork] = -Akkp1 / d;
               work[k+1+invd*ldwork] = -Akkp1 / d;
               k+=2;
            }
         }

         cut = n;
         while(cut > 0)
         {
            nnb = nb;
            if( cut <= nnb)
               nnb = cut;
            else
            {
               count = 0;
               for(i=cut-nnb;i<cut;i++)
                  if (bsdv[i] == 1)
                     count++;
               if(count%2 == 1)
                  nnb++;
            }

            cut-=nnb;

            for(i=0;i<cut;i++)
               for(j=0;j<nnb;j++)
                  work[i+j*ldwork] = A[i+(cut+j)*ldA];

            for(i=0;i<nnb;i++)
            {
               work[u11+i+1+i*ldwork] = one;
               for(j=0;j<i;j++)
                  work[u11+i+1+j*ldwork] = zero;
               for(j=i+1;j<nnb;j++)
                  work[u11+i+1+j*ldwork] = A[cut+i+(cut+j)*ldA];
            }

            i = 0;
            while(i<cut)
            {
               if(bsdv[i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[i+j*ldwork] = work[i+invd*ldwork]*work[i+j*ldwork];
                  i++;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u01_i_j = work[i+j*ldwork];
                     u01_ip1_j = work[i+1+j*ldwork];
                     work[i+j*ldwork] = work[i+invd*ldwork] * u01_i_j + work[i+(invd+1)*ldwork] * u01_ip1_j;
                     work[i+1+j*ldwork] = work[i+1+invd*ldwork] * u01_i_j + work[i+1+(invd+1)*ldwork] * u01_ip1_j;
                  }
                  i+=2;
               }
            }

            i = 0;
            while(i<nnb)
            {
               if(bsdv[cut+i] == 0)
               {
                  for(j=i;j<nnb;j++)
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork];
                  i++;
               }
               else
               {
                  for(j=i;j<nnb;j++)
                  {
                     u11_i_j = work[u11+i+1+j*ldwork];
                     u11_ip1_j = work[u11+i+2+j*ldwork];
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork] + work[cut+i+(invd+1)*ldwork] * work[u11+i+2+j*ldwork];
                     work[u11+i+2+j*ldwork] = work[cut+i+1+invd*ldwork] * u11_i_j + work[cut+i+1+(invd+1)*ldwork] * u11_ip1_j;
                  }
                  i+=2;
               }
            }

            TRMM<real_t>( 'L', 'U', 'T', 'U', nnb, nnb, one, A+cut+cut*ldA, ldA, work+u11+1, n+nb+1 );
            for(i=0;i<nnb;i++)
               for(j=i;j<nnb;j++)
                  A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];

            GEMM<real_t>( 'T', 'N', nnb, nnb, cut, one, A+cut*ldA, ldA, work, n+nb+1, zero, work+u11+1, n+nb+1 );
            for(i=0;i<nnb;i++)
               for(j=i;j<nnb;j++)
                  A[cut+i+(cut+j)*ldA] = A[cut+i+(cut+j)*ldA] + work[u11+i+1+j*ldwork];

            TRMM<real_t>( 'L', 'U', 'T', 'U', cut, nnb, one, A, ldA, work, n+nb+1 );

            for(i=0;i<cut;i++)
               for(j=0;j<nnb;j++)
                  A[i+(cut+j)*ldA] = work[i+j*ldwork];
         }

         i = 0;
         while(i<n)
         {
            ip = ipiv[i];
            if(bsdv[i] == 0)
            {
               if(i<ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i>ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
            }
            else
            {
               i++;
               if((i-1) < ip)
                  SYSWAPR( uplo, n, A, ldA, i-1, ip );
               if((i-1) > ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i-1 );
            }
            i++;
         }
      }
      else
      {

         TRTRI<real_t>( 'L', 'U', n, A, ldA, nb );

         k = n-1;
         while(k>=0)
         {
            if(bsdv[k] == 0)
            {
               work[k+invd*ldwork] = one / A[k+k*ldA];
               work[k+(invd+1)*ldwork] = zero;
               k--;
            }
            else
            {
               t = work[k-1];
               Ak = A[k-1+(k-1)*ldA] / t;
               Akp1 = A[k+k*ldA] / t;
               Akkp1 = work[k-1] / t;
               d = t * ( Ak * Akp1 - one );
               work[k-1+invd*ldwork] = Akp1 / d;
               work[k+invd*ldwork] = Ak / d;
               work[k+(invd+1)*ldwork] = -Akkp1 / d;
               work[k-1+(invd+1)*ldwork] = -Akkp1 / d;
               k-=2;
            }
         }

         cut = 0;
         while(cut < n)
         {
            nnb = nb;
            if(cut+nnb > n)
               nnb = n - cut;
            else
            {
               count = 0;
               for(i=cut;i<cut+nnb;i++)
                  if (bsdv[i] == 1)
                     count++;
               if(count%2 == 1)
                  nnb++;
            }

            for(i=0;i<n-cut-nnb;i++)
               for(j=0;j<nnb;j++)
                  work[i+j*ldwork] = A[cut+nnb+i+(cut+j)*ldA];

            for(i=0;i<nnb;i++)
            {
               work[u11+i+1+i*ldwork] = one;
               for(j=i+1;j<nnb;j++)
                  work[u11+i+1+j*ldwork] = zero;
               for(j=0;j<=i-1;j++)
                  work[u11+i+1+j*ldwork] = A[cut+i+(cut+j)*ldA];
            }

            i = n - 1 - cut - nnb;
            while(i>=0)
            {
               if(bsdv[cut+nnb+i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[i+j*ldwork] = work[cut+nnb+i+invd*ldwork]*work[i+j*ldwork];
                  i--;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u01_i_j = work[i+j*ldwork];
                     u01_ip1_j = work[i-1+j*ldwork];
                     work[i+j*ldwork] = work[cut+nnb+i+invd*ldwork] * u01_i_j + work[cut+nnb+i+(invd+1)*ldwork] * u01_ip1_j;
                     work[i-1+j*ldwork] = work[cut+nnb+i-1+(invd+1)*ldwork] * u01_i_j + work[cut+nnb+i-1+invd*ldwork] * u01_ip1_j;
                  }
                  i-=2;
               }
            }

            i = nnb - 1;
            while(i>=0)
            {
               if(bsdv[cut+i] == 0)
               {
                  for(j=0;j<nnb;j++)
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork];
                  i--;
               }
               else
               {
                  for(j=0;j<nnb;j++)
                  {
                     u11_i_j = work[u11+i+1+j*ldwork];
                     u11_ip1_j = work[u11+i+j*ldwork];
                     work[u11+i+1+j*ldwork] = work[cut+i+invd*ldwork] * work[u11+i+1+j*ldwork] + work[cut+i+(invd+1)*ldwork] * u11_ip1_j;
                     work[u11+i+j*ldwork] = work[cut+i-1+(invd+1)*ldwork] * u11_i_j + work[cut+i-1+invd*ldwork] * u11_ip1_j;
                  }
                  i-=2;
               }
            }

            TRMM<real_t>( 'L', 'L', 'T', 'U', nnb, nnb, one, A+cut+cut*ldA, ldA, work+u11+1, n+nb+1 );

            for(i=0;i<nnb;i++)
               for(j=0;j<=i;j++)
                  A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];

            if((cut+nnb) < n)
            {
               GEMM<real_t>( 'T', 'N', nnb, nnb, n-nnb-cut, one, A+cut+nnb+cut*ldA, ldA, work, n+nb+1, zero, work+u11+1, n+nb+1 );

               for(i=0;i<nnb;i++)
                  for(j=0;j<=i;j++)
                     A[cut+i+(cut+j)*ldA] = A[cut+i+(cut+j)*ldA] + work[u11+i+1+j*ldwork];

               TRMM<real_t>( 'L', 'L', 'T', 'U', n-nnb-cut, nnb, one, A+cut+nnb+(cut+nnb)*ldA, ldA, work, n+nb+1 );

               for(i=0;i<n-cut-nnb;i++)
                  for(j=0;j<nnb;j++)
                     A[cut+nnb+i+(cut+j)*ldA] = work[i+j*ldwork];
            }
            else
            {
               for(i=0;i<nnb;i++)
                  for(j=0;j<i;j++)
                     A[cut+i+(cut+j)*ldA] = work[u11+i+1+j*ldwork];
            }

            cut += nnb;
         }

         i = n - 1;
         while(i>=0)
         {
            ip = ipiv[i];
            if(bsdv[i] == 0)
            {
               if(i<ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i>ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
            }
            else
            {
               if(i < ip)
                  SYSWAPR( uplo, n, A, ldA, i, ip );
               if(i > ip)
                  SYSWAPR( uplo, n, A, ldA, ip, i );
               i--;
            }
            i--;
         }
      }

      delete [] work;
      return 0;
   }

}

#endif

