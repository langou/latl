//
// gebd2.h
// Linear Algebra Template Library
//
// Created by Philippe Theveny on 12/9/14.
// Copyright (c) 2014 University of Colorado Denver. All rights reserved.
//

#ifndef _gebd2_h
#define _gebd2_h

/// @file gebd2.h Reduces a general matrix to bidiagonal form using an unblocked algorithm.

#include "larf.h"
#include "larfg.h"

namespace LATL
{
  /// @brief Reduces a real general matrix to bidiagonal form using
  /// an unblocked algorithm.
  ///
  /// The real version of gebd2 reduces a real general m by n matrix A
  /// to upper or lower bidiagonal form B by an orthogonal tranformation:
  /// Q**T * A * P = B
  ///
  /// If m >= n, then B is upper diagonal, else B is lower diagonal.
  /// The matrices Q and P are represented as products of elementary
  /// reflectors:
  /// If m >= n,
  ///
  ///   Q = H(1) H(2) ... H(n) and P = G(1) G(2) ... G(n-1)
  ///
  /// Each H(i) and G(i) has the form:
  ///
  ///   H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T
  ///
  /// where tauq and taup are real scalars, and v and u are real vectors;
  /// v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
  /// u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
  /// tauq is stored in tauq(i) and taup in taup(i).
  ///
  /// If m < n,
  /// 
  ///   Q = H(1) H(2) ... H(m-1) and P = G(1) G(2) ... G(m)
  ///
  /// Each H(i) and G(i) has the form:
  ///
  ///   H(i) = I - tauq * v * v**T and G(i) = I - taup * u * u**T
  ///
  /// where tauq and taup are real scalars, and v and u are real vectors;
  /// v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
  /// u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
  /// tauq is stored in tauq(i) and taup in taup(i).
  ///
  /// @return 0 if success.
  /// @return -i if the ith argument is invalid.
  /// @tparam real_t Floating point type.
  /// @param[in] m The number of rows of the matrix A. m>=0
  /// @param[in] n The number of columns of the matrix A. n>=0
  /// @param[in,out] A pointer to real matrix A.
  ///  On entry, the m by n general matrix to be reduced.
  ///  On exit,
  ///  if m >= n, the diagonal and the first super diagonal are overwritten
  ///  with the upper bidiagonal matrix B; the elements below the diagonal,
  ///  with the array tauq, represent the orthogonal matrix Q as a product of
  ///  elementary reflectors, and the elements above the superdiagonal, with
  ///  the array taup, represent the orthogonal matrix P as a product of
  ///  elementary reflectors;
  ///  if m < n, the diagonal and the first subdiagonal are overwritten with
  ///  the lower bidiagonal matrix B; the elements below the first subdiagonal,
  ///  with the array tauq, represent the orthogonal matrix Q as a product of
  ///  elementary reflectors, and the elements above the diagonal, with the
  ///  array taup, represent the ortogonal matrix P as a product of elementary
  ///  reflectors.
  /// @param[in] ldA Column length of the matrix A. ldA>=max(1,m)
  /// @param[out] d The diagonal elements of the bidiagonal matrix B.
  ///  Array of dimension min(m,n).
  /// @param[out] e The off-diagonal elements of the bidiagonal matrix B.
  ///  Array of dimension min(m,n)-1.
  /// @param[out] tauq The scalar factors of the elementary reflectors which
  ///  represent the orthogonal matrix Q. Array of dimension min(m,n).
  /// @param[out] taup The scalar factors of the elementary reflectors which
  ///  represent the orthogonal matrix P. Array of dimension min(m,n).
  ///
  /// @ingroup COMP

  template<typename real_t>
  int GEBD2(int_t m, int_t n, real_t *A, int_t ldA, real_t *d, real_t *e,
	    real_t *tauq, real_t *taup)
  {
    int_t i;
    const real_t zero(0.0);
    const real_t one(1.0);
    using std::min;
    using std::max;
    
    if(m<0)
      return -1;
    else if(n<0)
      return -2;
    else if(ldA<m || ldA<1)
      return -4;

    real_t *work =new real_t[max(m,n)];

    if (m>=n)
    {
      // reduce to upper diagonal form
      for(i=0;i<n;i++)
      {
	LARFG<real_t>(m-i, A[i+i*ldA], &A[min(i+1,m-1)+i*ldA], 1, tauq[i]);
	d[i]=A[i+i*ldA];
	if(i<n-1)
	{
	  A[i+i*ldA]=one;
	  LARF<real_t>('L', m-i, n-i-1, &A[i+i*ldA], 1, tauq[i],
		       &A[i+(i+1)*ldA], ldA, work);
	  A[i+i*ldA]=d[i];
	  LARFG<real_t>(n-i-1, A[i+(i+1)*ldA],
			&A[i+min(i+2,n-1)*ldA], ldA, taup[i]);
	  e[i]=A[i+(i+1)*ldA];
	  A[i+(i+1)*ldA]=one;
	  LARF<real_t>('R', m-i-1, n-i-1, &A[i+(i+1)*ldA], ldA, taup[i],
		       &A[i+1+(i+1)*ldA], ldA, work);
	  A[i+(i+1)*ldA]=e[i];
	}
	else
	  taup[i]=zero;
      }
    }
    else
    {
      // reduce to lower diagonal form
      for(i=0;i<m;i++)
      {
	LARFG<real_t>(n-i, A[i+i*ldA], &A[i+min(i+1,n-1)*ldA], ldA, taup[i]);
	d[i]=A[i+i*ldA];
	if(i<m-1)
	{
	  A[i+i*ldA]=one;
	  LARF<real_t>('R', m-i-1, n-i, &A[i+i*ldA], ldA, taup[i],
		       &A[i+1+i*ldA], ldA, work);
	  A[i+i*ldA]=d[i];
	  LARFG<real_t>(m-i-1, A[i+1+i*ldA],
			&A[min(i+2,m-1)+i*ldA], 1, tauq[i]);
	  e[i]=A[i+1+i*ldA];
	  A[i+1+i*ldA]=one;
	  LARF<real_t>('L', m-i-1, n-i-1, &A[i+1+i*ldA], 1, tauq[i],
		       &A[i+1+(i+1)*ldA], ldA, work);
	  A[i+1+i*ldA]=e[i];
	}
	else
	  tauq[i]=zero;
      }
    }

    delete [] work;
    return 0;
  }
}

#endif /* _gebd2_h */
