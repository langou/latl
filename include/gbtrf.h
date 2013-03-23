//
//  gbtrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/5/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gbtrf_h
#define _gbtrf_h

/// @file gbtrf.h Computes an LU factorization of an m-by-n band matrix A using partial pivoting with row interchanges.

#include "imax.h"
#include "swap.h"
#include "imax.h"
#include "scal.h"
#include "copy.h"
#include "gemm.h"
#include "ger.h"
#include "laswp.h"
#include "trsm.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes an LU factorization of a real m-by-n band matrix A using partial pivoting.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid
   /// @return i+1 if the ith column of the matrix A has a zero pivot.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL of the array need not be set.  The jth column of A is stored in the jth column of the array AB as:  A(i, j) = AB(kL+kU+i-j, j)
   ///
   /// On exit, U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL+kU.  Multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of the array AB.  ldAB >= 2*kL+kU
   /// @param IPIV Integer array deminsion min(m, n).  The pivot indices, 0 <= i <= min(m, n), row i was interchanged with row IPIV(i).

   template< typename real_t>
   int_t GBTRF(const int_t m, const int_t n, const int_t kL, const int_t kU, real_t * const AB, const int_t ldAB, int_t * const IPIV)
   {
      if (m < 0)
         return -1;
      if (n < 0)
         return -2;
      if (kL < 0)
         return -3;
      if (kU < 0)
         return -4;
      if (ldAB < (2*kL+kU+1))
         return -6;

      if (m == 0 || n == 0)
         return 0;

      int_t kV = kL + kU, jU = 0, km, jp, info = 0;
      real_t * ABj = AB;
      const real_t zero(0.0);
      const real_t one(1.0);

      for (int_t j = kU+1; j < std::min(kV, n); ++j)
      {
         ABj = AB + ldAB*j;
         for (int_t i = kV-j+1; i < kL; ++i)
         {
            ABj[i] = zero;
         }
      }

      for (int_t j = 0; j < std::min(m,n); ++j)
      {
         if (j+kV < n)
         {
            ABj = AB + ldAB*(j+kV);
            for (int_t i = 0; i < kL; ++i)
            {
               ABj[i] = zero;
            }
         }

         ABj = AB+ldAB*j;
         km = std::min(kL, m-j-1);
         jp = LATL::IMAX(km+1, ABj+kV, 1);
         IPIV[j] = jp+j;
         if (ABj[kV+jp] != zero)
         {
            jU = std::max(jU, std::min(j+kU+jp, n-1));
            if (jp != 0)
            {
               LATL::SWAP(jU-j+1, ABj+kV+jp, ldAB-1, ABj+kV, ldAB-1);
            }
            if (km > 0)
            {
               LATL::SCAL(km, one/ABj[kV], ABj+kV+1, 1);
               if (jU > j)
                  LATL::GER(km, jU-j, -one, ABj+kV+1, 1, ABj+ldAB+kV-1, ldAB-1, ABj+ldAB+kV, ldAB-1);
            }
         }
         else
         {
            if (info == 0)
               info = j+1;
         }
      }

      return info;
   }

   /// @brief Computes an LU factorization of a complex m-by-n band matrix A using partial pivoting.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid
   /// @return i+1 if the ith column of the matrix A has a zero pivot.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL of the array need not be set.  The jth column of A is stored in the jth column of the array AB as:  A(i, j) = AB(kL+kU+i-j, j)
   ///
   /// On exit, U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL+kU.  Multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of the array AB.  ldAB >= 2*kL+kU
   /// @param IPIV Integer array deminsion min(m, n).  The pivot indices, 0 <= i <= min(m, n), row i was interchanged with row IPIV(i).

   template< typename real_t>
   int_t GBTRF(const int_t m, const int_t n, const int_t kL, const int_t kU, complex<real_t> * const AB, const int_t ldAB, int_t * const IPIV)
   {
      return LATL::GBTRF< complex<real_t> > (m, n, kL, kU, AB, ldAB, IPIV);
   }

   /// @brief Computes an LU factorization of a real m-by-n band matrix A using partial pivoting with row interchanges.
   ///
   /// This is the blocked version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kL The number of subdiagonals within the band of A.  kL >= 0.
   /// @param kU The number of superdiagonals within the band of A.  kU >= 0.
   /// @param AB Real matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of matrix A.  ldAB >= 2*kL+kU+1
   /// @param IPIV Permutation matrix size min(m,n).  On exit, row k of A was exchanged with IPIV[k].
   /// @param nb Block size.

   template< typename real_t>
   int_t GBTRF(const int_t m, const int_t n, const int_t kL, const int_t kU, real_t * const AB, const int_t ldAB, int_t * const IPIV, const int_t nb)
   {
      if (m < 0)
         return -1;
      if (n < 0)
         return -2;
      if (kL < 0)
         return -3;
      if (kU < 0)
         return -4;
      if (ldAB < (2*kL+kU+1))
         return -6;
      
      if (m == 0 || n == 0)
         return 0;
      
      int_t info = 0;
      if (nb > kL || nb <= 1)
         info = LATL::GBTRF(m, n, kL, kU, AB, ldAB, IPIV);
      else
      {
         const real_t zero(0.0);
         const real_t one(1.0);
         const int_t kV = kU+kL;
         int_t ldWork = nb+1;
         real_t * const Work13 = new real_t[(ldWork*nb)];
         real_t * const Work31 = new real_t[(ldWork*nb)];
         
         real_t * temp13 = Work13;
         for (int_t j = 0; j < nb; ++j)
         {
            for (int_t i = 0; i < j-1; ++i)
            {
               temp13[i] = zero;
            }
            temp13 += ldWork;
         }
         real_t * temp31 = Work31;
         for (int_t j =  0; j < nb; ++j)
         {
            for (int_t i = j+1; i < nb; ++i)
            {
               temp31[i] = zero;
            }
            temp31 += ldWork;
         }
         real_t * ABj = AB;
         for (int_t j = kU+1; j < std::min(kV, n); ++j)
         {
            ABj = AB+ldAB*j;
            for (int_t i = kV-j; i < kL; ++i)
            {
               ABj[i] = zero;
            }
         }
         
         int_t jU = 0, jb, i2, i3, km, nw, jp, jm;
         real_t * ABjj = AB;
         for (int_t j = 0; j < std::min(m, n); j+=nb)
         {
            ABj = AB+ldAB*j;
            jb = std::min(nb, std::min(m, n)-j);
            i2 = std::min(kL-jb, m-j-jb);
            i3 = std::min(jb, m-j-kL);
            for (int_t jj = j; jj < j+jb; ++jj)
            {
               if (jj+kV < n)
               {
                  ABjj = AB+ldAB*(jj+kV);
                  for (int_t i = 0; i < kL; ++i)
                  {
                     ABjj[i] = zero;
                  }
               }
               
               ABjj = AB+ldAB*jj;
               km = std::min(kL, m-jj-1);
               jp = LATL::IMAX(km+1, ABjj+kV, 1);
               IPIV[jj] = jp+jj-j;
               if (ABjj[kV+jp] != zero)
               {
                  jU = std::max(jU, std::min(jj+kU+jp, n-1));
                  if (jp != 0)
                  {
                     if ((jp+jj) < (j+kL))
                     {
                        LATL::SWAP(jb, ABj+kV+jj-j, ldAB-1, ABj+kV+jp+jj-j, ldAB-1);
                     }
                     else
                     {
                        LATL::SWAP(jj-j, ABj+kV+jj-j, ldAB-1, Work31+jp+jj-j-kL, ldWork);
                        LATL::SWAP(j+jb-jj, ABjj+kV, ldAB-1, ABjj+kV+jp, ldAB-1);
                     }
                  }
                  LATL::SCAL(km, one/ABjj[kV], ABjj+kV+1, 1);
                  
                  jm = std::min(jU, j+jb-1);
                  if (jm > jj)
                  {
                     LATL::GER(km, jm-jj, -one, ABjj+kV+1, 1, ABjj+ldAB+kV-1, ldAB-1, ABjj+ldAB+kV, ldAB-1);
                  }
               }
               else
               {
                  if (info == 0)
                     info = jj+1;
               }
               nw = std::min(jj-j+1, i3);
               if (nw > 0)
               {
                  LATL::COPY(nw, ABjj+kV+kL-jj+j, 1, Work31+ldWork*(jj-j), 1);
               }
            }
            
            if ((j+jb) < n)
            {
               int_t j2, j3, k2;
               
               j2 = std::min(jU-j+1, kV)-jb;
               j3 = std::max(0, jU-j-kV+1);
               LATL::laswp(j2, ABj+ldAB*jb+kV-jb, ldAB-1, 0, jb-1, IPIV+j);
               
               for (int_t i = j; i < j+jb; ++i)
               {
                  IPIV[i] += j;
               }
               
               k2 = j+jb+j2;  //k2 = fortran k2
               int_t ip, temp1 = k2;
               real_t temp;
               real_t * ABt = AB + ldAB*(k2);
               for (int_t i = 0; i < j3; ++i)
               {
                  //temp1 = fortran jj - 1
                  for (int_t ii = j+i; ii < j+jb; ++ii)
                  {
                     ip = IPIV[ii];
                     if (ip != ii)
                     {
                        temp = ABt[kV+ii-temp1];
                        ABt[kV+ii-temp1] = ABt[kV+ip-temp1];
                        ABt[kV+ip-temp1] = temp;
                     }
                  }
                  ABt += ldAB;
                  ++temp1;
               }
               
               if (j2 > 0)
               {
                  real_t * ABjjb = ABj + ldAB*jb;
                  LATL::TRSM('L', 'L', 'N', 'U', jb, j2, one, ABj+kV, ldAB-1, ABjjb+kV-jb, ldAB-1);
                  if (i2 > 0)
                     LATL::GEMM('N', 'N', i2, j2, jb, -one, ABj+kV+jb, ldAB-1, ABjjb+kV-jb, ldAB-1, one, ABjjb+kV, ldAB-1);
                  if (i3 > 0)
                     LATL::GEMM('N', 'N', i3, j2, jb, -one, Work31, ldWork, ABjjb+kV-jb, ldAB-1, one, ABjjb+kV+kL-jb, ldAB-1);
               }
               
               std::cout << "j3: " << j3 << std::endl;
               if (j3 > 0)
               {
                  temp13 = Work13;
                  real_t * ABjkv = ABj+ldAB*kV;
                  for (int_t jj = 0; jj < j3; ++jj)
                  {
                     for (int_t ii = jj; ii < jb; ++ii)
                     {
                        temp13[ii] = ABjkv[ii-jj];
                     }
                     temp13 += ldWork;
                     ABjkv += ldAB;
                  }
                  
                  ABjkv = ABj+ldAB*kV;
                  LATL::TRSM('L', 'L', 'N', 'U', jb, j3, one, ABj+kV, ldAB-1, Work13, ldWork);
                  
                  if (i2 > 0)
                  {
                     LATL::GEMM('N', 'N', i2, j3, jb, -one, ABj+kV+jb, ldAB-1, Work13, ldWork, one, ABjkv+jb+1, ldAB-1);
                  }
                  if (i3 > 0)
                  {
                  latl:GEMM('N', 'N', i3, j3, jb, -one, Work31, ldWork, Work13, ldWork, one, ABjkv+kL, ldAB-1);
                  }
                  
                  temp13 = Work13;
                  for (int_t jj = 0; jj < j3; ++jj)
                  {
                     for (int_t ii = jj; ii < jb; ++ii )
                     {
                        ABjkv[ii-jj] = temp13[ii];
                     }
                     ABjkv += ldAB;
                     temp13 += ldWork;
                  }
               }
            }
            else
            {
               for (int_t i = j; i < j+jb; ++i)
               {
                  IPIV[i] += j;
               }
            }
            
            for (int_t jj = j+jb-1; jj >= j; --jj)
            {
               jp = IPIV[jj] - jj;
               if (jp != 0)
               {
                  if (jp+jj < j+kL)
                  {
                     LATL::SWAP(jj-j, ABj+kV+jj-j, ldAB-1, ABj+kV+jp+jj-j, ldAB-1);
                  }
                  else
                  {
                     LATL::SWAP(jj-j, ABj+kV+jj-j, ldAB-1, Work31+jp+jj-j-kL, ldWork);
                  }
               }
               nw = std::min(i3, jj-j+1);
               if (nw > 0)
               {
                  LATL::COPY(nw, Work31+(jj-j)*ldWork, 1, AB+ldAB*jj+kV+kL-jj+j, 1);
               }
            }
         }
         
         delete [] Work13;
         delete [] Work31;
      }
      return info;
   }
   
   /// @brief Computes an LU factorization of a complex m-by-n band matrix A using partial pivoting with row interchanges.
   ///
   /// This is the blocked version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kL The number of subdiagonals within the band of A.  kL >= 0.
   /// @param kU The number of superdiagonals within the band of A.  kU >= 0.
   /// @param AB Complex matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of matrix A.  ldAB >= 2*kL+kU+1
   /// @param IPIV Permutation matrix size min(m,n).  On exit, row k of A was exchanged with IPIV[k].
   /// @param nb Block size.
   
   template< typename real_t>
   int_t GBTRF(const int_t m, const int_t n, const int_t kL, const int_t kU, complex<real_t> * const AB, const int_t ldAB, int_t * const IPIV, const int_t nb)
   {
      return LATL::GBTRF< complex<real_t> > (m, n, kL, kU, AB, ldAB, IPIV, nb);
   }
}

#endif
