//
//  gbtf2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/5/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gbtf2_h
#define _gbtf2_h

/// @file gbtf2.h Computes an LU factorization of an m-by-n band matrix A using partial pivoting with row interchanges.

#include "imax.h"
#include "swap.h"
#include "scal.h"
#include "ger.h"
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
   /// @param AB Real array size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.
   /// Rows 0 to kL of the array need not be set.  The jth column of A is stored in the jth column of the array AB as:
   ///
   ///           A(i, j) = AB(kL+kU+i-j, j)
   /// On exit, U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL+kU.
   /// Multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of the array AB.  ldAB >= 2*kL+kU
   /// @param ipiv Integer array deminsion min(m, n).  The pivot indices, 0 <= i <= min(m, n), row i was interchanged with
   /// row ipiv(i).
   /// @ingroup COMP

   template< typename real_t>
   int_t GBTF2(const int_t m, const int_t n, const int_t kL, const int_t kU, real_t * const AB, const int_t ldAB, int_t * const ipiv)
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

      using std::max;
      using std::min;
      int_t kV = kL + kU, jU = 0, km, jp, info = 0;
      real_t * ABj = AB;
      const real_t zero(0.0);
      const real_t one(1.0);

      for (int_t j = kU+1; j < min(kV, n); ++j)
      {
         ABj = AB + ldAB*j;
         for (int_t i = kV-j+1; i < kL; ++i)
         {
            ABj[i] = zero;
         }
      }

      for (int_t j = 0; j < min(m,n); ++j)
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
         km = min(kL, m-j-1);
         jp = LATL::IMAX(km+1, ABj+kV, 1);
         ipiv[j] = jp+j;
         if (ABj[kV+jp] != zero)
         {
            jU = max(jU, min(j+kU+jp, n-1));
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
   /// @param AB Complex array size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.
   /// Rows 0 to kL of the array need not be set.  The jth column of A is stored in the jth column of the array AB as:
   ///
   ///           A(i, j) = AB(kL+kU+i-j, j)
   /// On exit, U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL+kU.
   /// Multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB Column length of the array AB.  ldAB >= 2*kL+kU
   /// @param ipiv Integer array deminsion min(m, n).  The pivot indices, 0 <= i <= min(m, n), row i was interchanged
   /// with row ipiv(i).
   /// @ingroup COMP

   template< typename real_t>
   int_t GBTF2(const int_t m, const int_t n, const int_t kL, const int_t kU, complex<real_t> * const AB, const int_t ldAB, int_t * const ipiv)
   {
      return LATL::GBTF2< complex<real_t> > (m, n, kL, kU, AB, ldAB, ipiv);
   }

}

#endif
