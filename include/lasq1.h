//
//  lasq1.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq1_h
#define _lasq1_h

/// @file lasq1.h Computes the singular values of a square bidiagonal matrix.

#include <algorithm>
#include <cmath>
#include "copy.h"
#include "las2.h"
#include "lascl.h"
#include "lasq2.h"
#include "lasrt.h"
#include "lamch.h"
#include "latl.h"

namespace LATL
{
  /// @brief Computes the singular values of a real n-by-n bidiagonal matrix
  ///   with diagonal d and off-diagonal e. The singular values are computed
  ///   to high relative accuracy, in the absence of denormalization,
  ///   underflow, and overflow.
  ///
  /// See "Accurate singular values and differential qd algorithms" by K. V.
  /// Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
  /// 1994.
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] n     The number of rows and columns in the matrix. n>=0.
  /// @param[in,out] d Vector of dimension n.
  ///   On entry, d contains the diagonal elements of the bidiagonal matrix
  ///   whose SVD is desired.
  ///   On normal exit, d contains the singular values in decreasing order.
  /// @param[in,out] e Vector of dimension n.
  ///   On entry, e[1:n-1] contain the off-diagonal elements of the
  ///   bidiagonal matrix whose SVD is desired.
  ///   On exit, e is overwritten.
  /// @return 0 if success.
  /// @return -i if the ith argument had an illegal value.
  /// @return 1 if a split was marked by a positive value in e
  /// @return 2 if current block of z is not diagonalized after 100*n
  ///   iterations in inner while loop. On exit, d and e represent a matrix
  ///   which the calling subroutine could use to finish the computation, or
  ///   even feed back into LASQ1.
  /// @return 3 if the termination criterion of outer while loop is not met
  ///   (the program created more than n unreduced blocks).
  /// @ingroup AUX

  template<typename real_t>
  int_t LASQ1(int_t n, real_t *d, real_t *e)
  {
    using std::abs;
    using std::max;

    const real_t zero(0.0);
    const real_t eps=LAMCH<real_t>('P');
    const real_t safemin=LAMCH<real_t>('S');
    const real_t scale=sqrt(eps/safemin);

    int_t info=0;

    if(n<0)
      return -1;

    // Early return, when possible
    if(n==0)
      return 0;
    if(n==1)
    {
      d[0]=abs(d[1]);
      return 0;
    }

    real_t sigmax;
    real_t sigmin;
    if(n==2)
    {
      LAS2<real_t>(d[0],e[0],d[1],sigmin,sigmax);
      d[1]=sigmax;
      d[2]=sigmin;
      return 0;
    }

    // Estimate the largest singular value
    sigmax=zero;
    for(int_t i=0; i<n-1; i++)
    {
      d[i]=abs(d[i]);
      sigmax=max(sigmax, abs(e[i]));
    }
    d[n-1]=abs(d[n-1]);
    if(sigmax==zero)
    {
      // Early return, the matrix is already diagonal
      LASRT<real_t>('D',n,d);
      return 0;
    }
    for(int_t i=0; i<n; i++)
      sigmax=max(sigmax, d[i]);

    // Copy d and e into workspace in the Z format and scale
    real_t * work = new real_t[4*n];
    COPY(n,d,1,work,2);
    COPY(n-1,e,1,work+1,2);
    LASCL<real_t>('G',0,0,sigmax,scale,2*n-1,1,work,2*n-1);

    // Compute the q's and e's
    for(int_t i=0;i<2*n-1;i++)
      work[i]=work[i]*work[i];
    work[2*n-1]=0;

    info=LASQ2<real_t>(n,work);
    if(info==0)
    {
      for(int_t i=0;i<n;i++)
        d[i]=sqrt(work[i]);
      LASCL<real_t>('G',0,0,scale,sigmax,n,1,d,n);
    }
    else if(info==2)
    {
      // Maximum number of iterations exceeded
      // Move data from work into d and e so that the calling subroutine can
      // try to finish
      for(int_t i=0;i<n;i++)
      {
        d[i]=sqrt(work[2*i]);
        e[i]=sqrt(work[2*i+1]);
      }
      LASCL<real_t>('G',0,0,scale,sigmax,n,1,d,n);
      LASCL<real_t>('G',0,0,scale,sigmax,n,1,e,n);
    }

    delete [] work;
    return info;
  }
}

#endif
