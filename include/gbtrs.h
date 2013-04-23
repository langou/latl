//
//  gbtrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 4/17/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gbtrs_h
#define _gbtrs_h

/// @file gbtrs.h Solves a system of linear equations A * X = B.

#include "latl.h"
#include "gemv.h"
#include "ger.h"
#include "geru.h"
#include "swap.h"
#include "tbsv.h"
#include "lacgv.h"

namespace LATL
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///      A * X = B      or    A' * X = B
   ///
   /// where A is an n-by-n matrix stored in banded form and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid
   /// @tparam real_t Floating point type
   /// @param trans Specifies the form of the system of equations.
   ///      if trans == 'N', no transpose: A * X = B
   ///      if trans == 'T' or 'C', transpose: A' * X = B
   /// @param n The order of the matrix A.  n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.
   /// @param AB Real matrix.  On input, should contain the details of the LU factorization of the band matrix A as computed by GBTRF.  U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL + kU, and the multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB  Column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv The pivot indices from the LU factorization of A.
   /// @param B Real matrix size n-by-nrhs.  On entry, the right hand matrix B.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   
   template< typename real_t>
   int_t GBTRS(const char trans, const int_t n, const int_t kL, const int_t kU, const int_t nrhs, real_t * const AB, const int_t ldAB, int_t * const ipiv, real_t * const B, const int_t ldB)
   {
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -1;
      if ( n < 0)
         return -2;
      if (kL < 0)
         return -3;
      if (kU < 0)
         return -4;
      if (nrhs < 0)
         return -5;
      if (ldAB < (2*kL+kU+1))
         return -7;
      if (ldB < n)
         return -10;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      
      using std::min;
      const real_t one(1.0);
      int_t kD = kU + kL;
      bool low = (kL > 0)?1:0;
      real_t * ABj = AB, *Bi = B;
      int_t lm, l;
      
      if (trans == 'N' || trans == 'n')
      {
         if (low)
         {
            for (int_t j = 0; j < n-1; ++j)
            {
               lm = min(kL, n-j-1);
               l = ipiv[j];
               if ( l != j)
               {
                  LATL::SWAP(nrhs, B+l, ldB, B+j, ldB);
               }
               LATL::GER(lm, nrhs, -one, ABj+kD+1, 1, B+j, ldB, B+j+1, ldB);
               ABj += ldAB;
            }
         }
         for (int_t i = 0; i < nrhs; ++i)
         {
            LATL::TBSV('U', 'N', 'N', n, kL+kU, AB, ldAB, Bi, 1);
            Bi += ldB;
         }
      }
      else
      {
         for (int_t i = 0; i < nrhs; ++i)
         {
            LATL::TBSV('U','T','N', n, kL+kU, AB, ldAB, Bi, 1);
            Bi += ldB;
         }
         if (low)
         {
            for (int_t j = n-2; j >= 0; --j)
            {
               lm = min(kL, n-j-1);
               LATL::GEMV('T', lm, nrhs, -one, B+j+1, ldB, ABj+kD+1, 1, one, B+j, ldB);
               l = ipiv[j];
               if (l != j)
               {
                  LATL::SWAP(nrhs, B+l, ldB, B+j, ldB);
               }
               ABj += ldAB;
            }
         }
      }
      return 0;
   }
   
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B      or    A' * X = B     or A.' * X = B
   ///
   /// where A is an n-by-n matrix stored in banded form and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid
   /// @tparam real_t Floating point type
   /// @param trans Specifies the form of the system of equations.
   ///      if trans == 'N', no transpose:   A * X = B
   ///      if trans == 'C', conjugate transpose:  A' * X = B
   ///      if trans == 'T', transpose:   A.' * X = B
   /// @param n The order of the matrix A.  n >= 0
   /// @param kL Number of subdiagonals within the band of A.  kL >= 0
   /// @param kU Number of superdiagonals within the band of A.  kU >= 0
   /// @param nrhs Number of columns of the matrix B.
   /// @param AB Complex matrix.  On input, should contain the details of the LU factorization of the band matrix A as computed by GBTRF.  U is stored as an upper triangular band matrix with kL + kU superdiagonals in rows 0 to kL + kU, and the multipliers used during the factorization are stored in rows kL+kU+1 to 2*kL+kU.
   /// @param ldAB  Column length of the matrix A.  ldAB >= 2*kL+kU+1
   /// @param ipiv The pivot indices from the LU factorization of A.
   /// @param B Complex matrix size n-by-nrhs.  On entry, the right hand matrix B.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   
   template< typename real_t>
   int_t GBTRS(const char trans, const int_t n, const int_t kL, const int_t kU, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, int_t * const ipiv, complex<real_t> * const B, const int_t ldB)
   {
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -1;
      if ( n < 0)
         return -2;
      if (kL < 0)
         return -3;
      if (kU < 0)
         return -4;
      if (nrhs < 0)
         return -5;
      if (ldAB < (2*kL+kU+1))
         return -7;
      if (ldB < n)
         return -10;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      
      using std::min;
      const complex<real_t> one(1.0);
      int_t kD = kU + kL;
      bool low = (kL > 0)?1:0;
      complex<real_t> * ABj = AB, *Bi = B;
      int_t lm, l;
      
      if (trans == 'N' || trans == 'n')
      {
         if (low)
         {
            for (int_t j = 0; j < n-1; ++j)
            {
               lm = min(kL, n-j-1);
               l = ipiv[j];
               if ( l != j)
               {
                  LATL::SWAP(nrhs, B+l, ldB, B+j, ldB);
               }
               LATL::GERU(lm, nrhs, -one, ABj+kD+1, 1, B+j, ldB, B+j+1, ldB);
               ABj += ldAB;
            }
         }
         for (int_t i = 0; i < nrhs; ++i)
         {
            LATL::TBSV('U', 'N', 'N', n, kD, AB, ldAB, Bi, 1);
            Bi += ldB;
         }
      }
      else if (trans == 'T' || trans == 't')
      {
         for (int_t i = 0; i < nrhs; ++i)
         {
            LATL::TBSV('U','T','N', n, kD, AB, ldAB, Bi, 1);
            Bi += ldB;
         }
         if (low)
         {
            for (int_t j = n-2; j >= 0; --j)
            {
               lm = min(kL, n-j-1);
               LATL::GEMV('T', lm, nrhs, -one, B+j+1, ldB, ABj+kD+1, 1, one, B+j, ldB);
               l = ipiv[j];
               if (l != j)
               {
                  LATL::SWAP(nrhs, B+l, ldB, B+j, ldB);
               }
               ABj += ldAB;
            }
         }
      }
      else
      {
         for (int_t i = 0; i < nrhs; ++i)
         {
            LATL::TBSV('U', 'C', 'N', n, kD, AB, ldAB, Bi, 1);
            Bi += ldB;
         }
         if (low)
         {
            for (int_t j = n-2; j >= 0; --j)
            {
               lm = min(kL, n-j-1);
               LATL::LACGV(nrhs, B+j, ldB);
               LATL::GEMV('C', lm, nrhs, -one, B+j+1, ldB, ABj+kD+1, 1, one, B+j, ldB);
               LATL::LACGV(nrhs, B+j, ldB);
               l = ipiv[j];
               if (l != j)
               {
                  LATL::SWAP(nrhs, B+l, ldB, B+j, ldB);
               }
               ABj += ldAB;
            }
         }
      }
      return 0;
   }
}


#endif
