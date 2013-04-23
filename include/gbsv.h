//
//  gbsv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 4/22/13.
//
//

#ifndef _gbsv_h
#define _gbsv_h

#include "latl.h"
#include "gbtrs.h"
#include "gbtrf.h"

/// @file gbsv.h Solves a system of equations A * X = B.

namespace LATL
{
   /// @brief Solves a real system of equations A * X = B where A is a real n-by-n matrix stored in banded form, and X and B are n-by-nrhs matrices.
   ///
   /// This version will use the unblocked algorithm GBTF2 to compute the LU factorization of A.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of the matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A. n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= n
   /// @param AB Real matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv Permutation matrix size n.
   /// @param B Real matrix size n-by-nrhs.  On exit, the solution matrix B.
   /// @param ldB Column length of the matrix B.
   
   template< typename real_t>
   int_t GBSV(const int_t n, const int_t kL, const int_t kU, const int_t nrhs, real_t * const AB, const int_t ldAB, int_t * const ipiv, real_t * const B, const int_t ldB)
   {
      if (n < 0)
         return -1;
      if (kL < 0)
         return -2;
      if (kU < 0)
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < 2*kL+kU+1)
         return -6;
      if (ldB < n)
         return -9;
      
      int_t info = LATL::GBTF2(n, n, kL, kU, AB, ldAB, ipiv);
      if (info == 0)
      {
         return LATL::GBTRS('N', n, kL, kU, nrhs, AB, ldAB, ipiv, B, ldB);
      }
   }
   
   /// @brief Solves a complex system of equations A * X = B where A is a complex n-by-n matrix stored in banded form, and X and B are n-by-nrhs matrices.
   ///
   /// This version will use the unblocked algorithm GBTF2 to compute the LU factorization of A.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of the matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A. n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= n
   /// @param AB Complex matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv Permutation matrix size n.
   /// @param B Complex matrix size n-by-nrhs.  On exit, the solution matrix B.
   /// @param ldB Column length of the matrix B.
   
   template< typename real_t>
   int_t GBSV(const int_t n, const int_t kL, const int_t kU, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, int_t * const ipiv, complex<real_t> * const B, const int_t ldB)
   {
      return LATL::GBSV< complex<real_t> >(n, kL, kU, nrhs, AB, ldAB, ipiv, B, ldB);
   }
   
   /// @brief Solves a real system of equations A * X = B where A is a real n-by-n matrix stored in banded form, and X and B are n-by-nrhs matrices.
   ///
   /// This version will use the blocked algorithm GBTRF to compute the LU factorization of A.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of the matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A. n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= n
   /// @param AB Real matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv Permutation matrix size n.
   /// @param B Real matrix size n-by-nrhs.  On exit, the solution matrix B.
   /// @param ldB Column length of the matrix B.
   /// @param nb Block size.
   
   template< typename real_t>
   int_t GBSV(const int_t n, const int_t kL, const int_t kU, const int_t nrhs, real_t * const AB, const int_t ldAB, int_t * const ipiv, real_t * const B, const int_t ldB, const int_t nb)
   {
      if (n < 0)
         return -1;
      if (kL < 0)
         return -2;
      if (kU < 0)
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < 2*kL+kU+1)
         return -6;
      if (ldB < n)
         return -9;
      
      int_t info = LATL::GBTRF(n, n, kL, kU, AB, ldAB, ipiv, nb);
      if (info == 0)
      {
         return LATL::GBTRS('N', n, kL, kU, nrhs, AB, ldAB, ipiv, B, ldB);
      }
   }
   
   /// @brief Solves a complex system of equations A * X = B where A is a complex n-by-n matrix stored in banded form, and X and B are n-by-nrhs matrices.
   ///
   /// This version will use the blocked algorithm GBTRF to compute the LU factorization of A.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of the matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A. n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= n
   /// @param AB Complex matrix size ldAB-by-n.  On entry, the matrix A in band storage, in rows kL to 2*kL+kU.  Rows 0 to kL need not be set.  On exit, the factor U is stored as an upper triangular band matrix with kL+kU superdiagonals in rows 0 to kL+kU, and the multipliers using during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv Permutation matrix size n.
   /// @param B Complex matrix size n-by-nrhs.  On exit, the solution matrix B.
   /// @param ldB Column length of the matrix B.
   /// @param nb Block size.
   
   template< typename real_t>
   int_t GBSV(const int_t n, const int_t kL, const int_t kU, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, int_t * const ipiv, complex<real_t> * const B, const int_t ldB, const int_t nb)
   {
      return LATL::GBSV< complex<real_t> >(n, kL, kU, nrhs, AB, ldAB, ipiv, B, ldB, nb);
   }
}

#endif
