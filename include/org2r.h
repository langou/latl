//
//  org2r.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 12/16/14.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _org2r_h
#define _org2r_h

/// @file org2r.h Generates all or part of the orthogonal matrix Q from a QR
/// factorization determined by geqrf (unblocked algorithm)

#include "larf.h"
#include "scal.h"
#include "latl.h"

namespace LATL
{
  /// @brief Generates a m-by-n real matrix Q with orthogonal columns.
  ///
  /// The matrix Q is represented as a product of elementary reflectors
  ///
  ///         Q = H_1 H_2 . . . H_k
  /// Each H_i has the form
  ///
  ///         H_i = I - tau * v * v'
  /// where tau is a scalar, and v is a vector with
  ///
  ///         v[0] through v[i-1] = 0; v[i]   = 1
  /// with v[i+1] through v[m-1] stored on exit below the diagonal
  /// in the ith column of A, and tau in tau[i].
  /// @tparam real_t Floating point type.
  /// @return 0 if success
  /// @return -i if the ith argument is invalid.
  /// @param[in] m The number of rows of the matrix A. m>=0
  /// @param[in] n The number of columns of the matrix A. n>=0
  /// @param]in] k The number of elementary reflectors whose product defines
  ///  the matrix Q. n>=k>=0
  /// @param[in,out] A Real m-by-n matrix.
  ///  On entry, the i-th column must contains the vector which defines the
  ///  elementary reflector H(i), for i=1,2,...,k, as returned by GEQRF in the
  ///  first k columns of its array argument A.
  ///  On exit, the m-by-n matrix Q.
  /// @param[in] ldA Column length of the matrix A.  ldA>=m
  /// @param[out] tau Real vector of length k.
  ///  The scalar factors of the elementary reflectors.
  /// @ingroup COMP

  template<typename real_t>
  int ORG2R(int_t m, int_t n, int_t k, real_t *A, int_t ldA, real_t *tau)
  {
    const real_t zero(0.0);
    const real_t one(1.0);
     
    if (m<0)
      return -1;
    else if (n<0 || n>m)
      return -2;
    else if (k<0 || k>n)
      return -3;
    else if (ldA<m)
      return -5;

    if (n<=0)
      // quick return
      return 0;

    real_t *work=new real_t[n];
     
    for (int_t j=k;j<n;j++)
    {
      for (int_t l=0;l<m;l++)
	A[l+j*ldA]=zero;

      A[j+j*ldA]=one;
    }

    for (int_t i=k-1; i>=0; i--)
    {
      if (i<n-1)
      {
	A[i+i*ldA]=one;
	LARF<real_t>('L', m-i, n-i-1, &A[i+i*ldA], 1, tau[i],
		     &A[i+(i+1)*ldA], ldA, work);
      }
      if (i<m-1)
	SCAL<real_t>(m-i-1, -tau[i], &A[i+1+i*ldA], 1);
      A[i+i*ldA]=one-tau[i];

      for (int_t l=0; l<i; l++)
	A[l+i*ldA]=zero;
    }
    
    delete [] work;
    return 0;
  }
}
#endif
