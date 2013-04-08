//
//  gebnd.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/7/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _gebnd_h
#define _gebnd_h

/// @file gebnd.h Copies a band matrix from general format to banded format.


#include "latl.h"

namespace LATL
{
   /// @brief Copies a band matrix from general format to banded format.
   ///
   /// The matrix bands are stored in banded format as rows as follows.
   /// The diagonal is stored as row ku, starting in column 0.  The first superdiagonal is stored starting in
   /// column 1 of row ku-1, and the second superdiagonal in column 2 of row ku-2, and so on.  The first subdiagonal
   /// is stored starting in column 0 of row ku+1, and the second subdiagonal in column 0 of row ku+2, and so on.
   /// As an example, consider the following banded 5-by-5 matrix with ku=2 and kl=1.  On the left is the matrix
   /// in standard form, and on the right the same matrix in banded form.
   ///
   ///          (general)            (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ p i e b . ]        [ . d e f g ]
   ///        [ . q j f c ]        [ h i j k l ]
   ///        [ . . r k g ]        [ p q r s . ]
   ///        [ . . . s l ]
   ///
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param n Order of band matrix. n>0.
   /// @param[in] A Real matrix of order n.  A is assumed to be a band matrix stored in general form.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param kl Lower band depth.  n>kl>=0.
   /// @param ku Upper band depth.  n>ku>=0, n>ku+kl.
   /// @param[out] B Real m-by-n matrix containing the bands of A, where m=1+kl+ku.
   /// Each row of B contains a band of A, starting with the super diagonal ku in the first row.
   /// @param ldB Column length of the matrix B.  ldB>=1+kl+ku.
   /// @ingroup AUX

   template<typename real_t>
   int GEBND(int_t n,real_t *A,int_t ldA,int_t kl,int_t ku,real_t *B,int_t ldB)
   {
      const real_t zero(0.0);
      if(n<0)
         return -1;
      else if(ldA<n)
         return -3;
      else if(kl>n-1)
         return -4;
      else if(ku+kl+1>n)
         return -5;
      else if(ldB<ku+kl+1)
         return -6;

      for(int_t j=0;j<n;j++)
         B[ku+j*ldB]=A[j+j*ldA];

      for(int_t i=ku;i>0;--i)
      {
         for(int_t j=0;j<i;j++)
            B[ku-i+j*ldB]=zero;
         for(int_t j=i;j<n;j++)
            B[ku-i+j*ldB]=A[j-i+j*ldA];
      }

      for(int_t i=1;i<=kl;i++)
      {
         for(int_t j=0;j<n-i;j++)
            B[ku+i+j*ldB]=A[i+j+j*ldA];
         for(int_t j=n-i;j<n;j++)
            B[ku+i+j*ldB]=zero;
      }

      return 0;
   }

   /// @brief Copies a band matrix from general format to banded format.
   ///
   /// The matrix bands are stored in banded format as rows as follows.
   /// The diagonal is stored as row ku, starting in column 0.  The first superdiagonal is stored starting in
   /// column 1 of row ku-1, and the second superdiagonal in column 2 of row ku-2, and so on.  The first subdiagonal
   /// is stored starting in column 0 of row ku+1, and the second subdiagonal in column 0 of row ku+2, and so on.
   /// As an example, consider the following banded 5-by-5 matrix with ku=2 and kl=1.  On the left is the matrix
   /// in standard form, and on the right the same matrix in banded form.
   ///
   ///          (general)            (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ p i e b . ]        [ . d e f g ]
   ///        [ . q j f c ]        [ h i j k l ]
   ///        [ . . r k g ]        [ p q r s . ]
   ///        [ . . . s l ]
   ///
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param n Order of band matrix. n>0.
   /// @param[in] A Complex matrix of order n.  A is assumed to be a band matrix stored in general form.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param kl Lower band depth.  n>kl>=0.
   /// @param ku Upper band depth.  n>ku>=0, n>ku+kl.
   /// @param[out] B Complex m-by-n matrix containing the bands of A, where m=1+kl+ku.
   /// Each row of B contains a band of A, starting with the super diagonal ku in the first row.
   /// @param ldB Column length of the matrix B.  ldB>=1+kl+ku.
   /// @ingroup AUX

   template<typename real_t>
   int GEBND(int_t n,complex<real_t> *A,int_t ldA,int_t kl,int_t ku,complex<real_t> *B,int_t ldB)
   {
      return GEBND< complex<real_t> >(n,A,ldA,kl,ku,B,ldB);
   }
}
#endif
