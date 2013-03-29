//
//  laqge.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/5/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _laqge_h
#define _laqge_h

/// @file laqge.h  Equilibrates a general m-by-n matrix A using the row and column scaling factors in the vectors R and C.

#include <limits>
#include "latl.h"

namespace LATL
{
   /// @brief Equilibrates a general m-by-n matrix A using the row and column scaling factors in the vectors R and C.
   ///
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in the matrix A.
   /// @param n Number of columns in the matrix A.
   /// @param A Real array size ldA-by-n.  On entry, contains the m-by-n matrix A.  On exit, the equilibrated matrix.
   /// @param ldA Column size of the matrix A.
   /// @param R Real array of size m.  On entry, should contain the row scale factors for A.
   /// @param C Real array of size n.  On entry, should contain the column scale factors for A.
   /// @param rowcnd Real number; the ratio of the smallest R[i] to the largest R[i].
   /// @param colcnd Real number; the ratio of the smallest C[i] to the largest C[i].
   /// @param amax Real number; absolute valu eof the largest matrix entry.
   /// @param equed Character.  On exit, indicates the form of equilibration that was done.
   ///
   ///     'N': No equilibration
   ///     'R': Row equilibration; A has been premultiplied by diag(R).
   ///     'C': Column equilibration; A has been postmultiplied by diag(C).
   ///     'B': Both row and column equilibration; A has been replaced by diag(R) * A * diag(C).
   /// @ingroup AUX
   
   template< typename real_t>
   void LAQGE(const int_t m, const int_t n, real_t * const A, const int_t ldA, real_t * const R, real_t * const C, const real_t rowcnd, const real_t colcnd, const real_t amax, char &equed)
   {
      if ( m <= 0 or n <= 0)
      {
         equed = 'N';
         return;
      }
      
      const real_t one(1.0);
      const real_t thresh(0.1);
      real_t small = numeric_limits<real_t>::min()/numeric_limits<real_t>::epsilon();
      real_t large = one/small;
      
      real_t * Aj = A;
      
      if (rowcnd >= thresh && amax >= small && amax <= large)
      {
         if (colcnd >= thresh)
         {
            equed = 'N';
            return;
         }
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Aj[i] *= C[j];
            }
            Aj += ldA;
         }
         equed = 'C';
         return;
      }
      if (colcnd >= thresh)
      {
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Aj[i] *= R[i];
            }
            Aj += ldA;
         }
         equed = 'R';
         return;
      }
      for (int_t j = 0; j < n; ++j)
      {
         for (int_t i = 0; i< m; ++i)
         {
            Aj[i] = C[j] * R[i] * Aj[i];
         }
         Aj += ldA;
      }
      equed = 'B';
      return;
   }
   
   /// @brief Equilibrates a general m-by-n matrix A using the row and column scaling factors in the vectors R and C.
   ///
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in the matrix A.
   /// @param n Number of columns in the matrix A.
   /// @param A Complex array size ldA-by-n.  On entry, contains the m-by-n matrix A.  On exit, the equilibrated matrix.
   /// @param ldA Column size of the matrix A.
   /// @param R Real array of size m.  On entry, should contain the row scale factors for A.
   /// @param C Real array of size n.  On entry, should contain the column scale factors for A.
   /// @param rowcnd Real number; the ratio of the smallest R[i] to the largest R[i].
   /// @param colcnd Real number; the ratio of the smallest C[i] to the largest C[i].
   /// @param amax Real number; absolute valu eof the largest matrix entry.
   /// @param equed Character.  On exit, indicates the form of equilibration that was done.
   ///
   ///     'N': No equilibration
   ///     'R': Row equilibration; A has been premultiplied by diag(R).
   ///     'C': Column equilibration; A has been postmultiplied by diag(C).
   ///     'B': Both row and column equilibration; A has been replaced by diag(R) * A * diag(C).
   /// @ingroup AUX

   template< typename real_t>
   void LAQGE(const int_t m, const int_t n, complex<real_t> * const A, const int_t ldA, real_t * const R, real_t * const C, const real_t rowcnd, const real_t colcnd, const real_t amax, char &equed)
   {
      
      if ( m <= 0 or n <= 0)
      {
         equed = 'N';
         return;
      }
      
      const real_t one(1.0);
      const real_t thresh(.1);
      real_t small = numeric_limits<real_t>::min()/numeric_limits<real_t>::epsilon();
      real_t large = one/small;
      
      complex<real_t> * Aj = A;
      
      if (rowcnd >= thresh && amax >= small && amax <= large)
      {
         if (colcnd >= thresh)
         {
            equed = 'N';
            return;
         }
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Aj[i] *= C[j];
            }
            Aj += ldA;
         }
         equed = 'C';
         return;
      }
      if (colcnd >= thresh)
      {
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Aj[i] *= R[i];
            }
            Aj += ldA;
         }
         equed = 'R';
         return;
      }
      for (int_t j = 0; j < n; ++j)
      {
         for (int_t i = 0; i< m; ++i)
         {
            Aj[i] = C[j] * R[i] * Aj[i];
         }
         Aj += ldA;
      }
      equed = 'B';
      return;
   }
   
}

#endif
