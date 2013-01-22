//
//  lag2x.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lag2x_h
#define _lag2x_h

/// @file lag2x.h Converts a real matrix to higher precision.

#include "latl.h"

namespace latl
{
   /// @brief Convert a real matrix to higher precision.
   ///
   /// The matrix A of is copied to the matrix B.
   ///
   ///        B := A
   /// @return 0 if success.
   /// @return -i if ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @tparam xreal_t Real floating point type of higher precision.
   /// @param m Number of rows of the matrix A.
   /// @param n Number of columns of the matrix A.
   /// @param A Pointer to matrix A.
   /// @param ldA Column length of matrix A.  ldA>=m
   /// @param B Pointer to matrix B.
   /// @param ldB Column length of matrix B.  ldB>=m
   /// @ingroup MAT
   
   template<typename real_t,typename xreal_t>
   int lag2x(int_t m,int_t n,real_t *A,int_t ldA,xreal_t *B,int_t ldB)
   {
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(ldB<m)
         return -6;

      for(int_t j=0;j<n;j++)
      {
         for(int_t i=0;i<m;i++)
            B[i]=A[i];
         A+=ldA;
         B+=ldB;
      }
      return info;
   }

   /// @brief Convert a complex matrix to higher precision.
   ///
   /// The matrix A of is copied to the matrix B.
   ///
   ///        B := A
   /// @return 0 if success.
   /// @return -i if ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @tparam xreal_t Real floating point type of higher precision.
   /// @param m Number of rows of the matrix A.
   /// @param n Number of columns of the matrix A.
   /// @param A Pointer to complex matrix A.
   /// @param ldA Column length of matrix A.  ldA>=m
   /// @param B Pointer to complex matrix B.
   /// @param ldB Column length of matrix B.  ldB>=m
   /// @ingroup MAT
   
   template<typename real_t,typename xreal_t>
   int lag2x(int_t m,int_t n,complex<real_t> *A,int_t ldA,complex<xreal_t> *B,int_t ldB)
   {
      return lag2x< complex<real_t>, complex<xreal_t> >(m,n,A,ldA,B,ldB);
   }
}

#endif
