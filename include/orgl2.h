//
//  orgl2.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 12/16/14.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _orgl2_h
#define _orgl2_h

/// @file orgl2.h Generates all or part of the orthogonal matrix Q from a LQ
/// factorization determined by gelqf

#include "larf.h"
#include "scal.h"
#include "latl.h"

namespace LATL
{
  /// @brief Generates a m-by-n real matrix Q with orthogonal rows.
  ///
  /// The matrix Q is represented as a product of elementary reflectors
  ///
  ///         Q = H_1 H_2 . . . H_k
  /// Each H_i has the form
  ///
  ///         H_i = I - tau * u * u'
  /// where tau is a scalar, and u is a vector with
  ///
  ///         u[0] through u[i-1] = 0; u[i]   = 1
  /// with u[i+1] through u[m-1] stored on exit above the diagonal
  /// in the ith row of A, and tau in tau[i].
  /// @tparam real_t Floating point type.
  /// @return 0 if success
  /// @return -i if the ith argument is invalid.
  /// @param[in] m The number of rows of the matrix Q. m>=0
  /// @param[in] n The number of columns of the matrix Q. n>=m
  /// @param]in] k The number of elementary reflectors whose product defines
  ///  the matrix Q. m>=k>=0
  /// @param[in,out] A Real ldA-by-n matrix.
  ///  On entry, the i-th row must contains the vector which defines the
  ///  elementary reflector H(i), for i=1,2,...,k, as returned by GELQF in the
  ///  first k rows of its array argument A.
  ///  On exit, the m-by-n matrix Q.
  /// @param[in] ldA Column length of the matrix A.  ldA>=m
  /// @param[out] tau Real vector of length k.
  ///  The scalar factors of the elementary reflectors.
  /// @ingroup COMP

  template<typename real_t>
  int ORGL2(int_t m, int_t n, int_t k, real_t *A, int_t ldA, real_t *tau)
  {
    const real_t zero(0.0);
    const real_t one(1.0);
     
    if (m<0)
      return -1;
    else if (n<m)
      return -2;
    else if (k<0 || m<k)
      return -3;
    else if (ldA<m)
      return -5;

    if (m<=0)
      // quick return
      return 0;

    real_t *work=new real_t[m];

    if (k<m)
      for (int_t j=0;j<n;j++)
      {
	for (int_t l=k;l<m;l++)
	  A[l+j*ldA]=zero;
	if (j>k-1 && j<m)
	  A[j+j*ldA]=one;
      }

    for (int_t i=k-1; i>=0; i--)
    {
      if (i<n-1)
      {
      	if (i<m-1)
      	{
	  A[i+i*ldA]=one;
	  LARF<real_t>('R', m-i-1, n-i, &A[i+i*ldA], ldA, tau[i],
		       &A[i+1+i*ldA], ldA, work);
	}
	SCAL<real_t>(n-i-1, -tau[i], &A[i+(i+1)*ldA], ldA);
      }
      A[i+i*ldA]=one-tau[i];

      for (int_t l=0; l<i; l++)
	A[i+l*ldA]=zero;
    }
    
    delete [] work;
    return 0;
  }
}
#endif
