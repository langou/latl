//
//  lasq3.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq3_h
#define _lasq3_h

/// @file lasq3.h Checks for deflation, computes a shift and calls dqds.

#include <algorithm>
#include <cmath>
#include "lasq4.h"
#include "lasq5.h"
#include "lasq6.h"
#include "lamch.h"
#include "latl.h"

namespace LATL
{
  /// @brief Checks for deflation, computes a shift (tau) and calls dqds.
  ///   In case of failure, it changes shifts, and tries again until output
  ///   is positive.
  ///
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] i0        First index.
  /// @param[in,out] n0    Last index.
  /// @param[in,out] z     Real array of length 4*n. z holds the qd array.
  /// @param[in,out] pp    0 for ping, 1 for pong. pp=2 indicates that
  ///   flipping was applied to the z array and that the initial tests for
  ///   deflation should not be performed.
  /// @param[out] dmin     Minimum value of d.
  /// @param[out] sigma    The sum of shifts used in the current segment.
  /// @param[in,out] desig Lower order part of sigma.
  /// @param[in] qmax      Maximum value of q.
  /// @param[out] nfail    Number of times the shif was too big.
  /// @param[out] iter     Number of iterations.
  /// @param[out] ndiv     Number of divisions.
  /// @param[in,out] ttype Shift type.
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] dmin1 Minimum value of d, excluding d[n0-1].
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] dmin2 Minimum value of d, excluding d[n0-1] and d[n0-2].
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] dn    d[n0-1], the last value of d.
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] dn1   d[n0-2].
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] dn2   d[n0-3].
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] g     fraction of dmin used.
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @param[in,out] tau   The shift.
  ///   Passed as an argument in order to save value between calls to LASQ3.
  /// @ingroup AUX
  template<typename real_t>
  void LASQ3(int_t i0, int_t &n0, real_t *z, int_t &pp, real_t &dmin,
             real_t &sigma, real_t &desig, int_t qmax,
             int_t &nfail, int_t &iter, int_t &ndiv, int_t &ttype,
             real_t &dmin1, real_t &dmin2, real_t &dn, real_t &dn1,
             real_t &dn2, real_t &g, real_t &tau)
  {
    using std::abs;
    using std::min;
    using std::max;
    using std::sqrt;

    if (n0<i0)
      return;

    const real_t zero(0.0);
    const real_t one(1.0);
    const real_t two(2.0);
    const real_t half(0.5);
    const real_t quarter(0.25);
    const real_t hundred(100.0);
    const real_t cbias(1.5);
    const real_t eps=LAMCH<real_t>('P');
    const real_t tol=eps*hundred;
    const real_t tol2=tol*tol;

    int_t j, n0in;
    real_t t;

    // z is in qd-format
    // z[4(i-1)]   is qi, z[4(i-1)+1] is qqi,
    // z[4(i-1)+2] is ei, z[4(i-1)+3] is eei

    // check for deflation
    while(n0>i0)
    {
      if(n0>i0+1)
      {
        if(z[4*(n0-2)+2+pp]<=tol2*(sigma+z[4*(n0-1)+pp])
           || z[4*(n0-2)+3-pp]<=tol2*z[4*(n0-1)+pp])
        {
          // e(n0-1) is negligible, 1 eigenvalue
          z[4*(n0-1)]=z[4*(n0-1)+pp]+sigma;
          no--;
          continue;
        }
        else if(z[4*(n0-3)+2+pp]>tol2*sigma
                && z[4*(n0-3)+3-pp]>tol2*z[4*(n0-3)+pp])
        {
        // neither e(n0-1) nor e(n0-2) is negligible, no deflation

        // warning: LAPACK 3.4.2 code and the article referenced above do not
        // match.  Here, we follow LAPACK 3.4.2 in the criterion for
        // neglecting e(n0-2).
        break;
        }
      }

      // either n0=i0+1, or e(n0-2) is negligible
      if(z[4*(n0-1)+pp]>z[4*(n0-2)+pp])
      {
        // put q's in decreasing order
        t=z[4*(n0-1)+pp];
        z[4*(n0-1)+pp]=z[4*(n0-2)+pp];
        z[4*(n0-2)+pp]=t;
      }

      // deflate for two eigenvalues following Rutishauser's algorithm
      if(z[4*(n0-2)+2+pp]>tol2*z[4*(n0-1)+pp])
      {
        t=half*((z[4*(n0-2)+pp]-z[4*(n0-1)+pp])+z[4*(n0-2)+2+pp]);
        if(t!=zero)
        {
          real_t s;
          s=z[4*(n0-1)+pp]*(z[4*(n0-2)+2+pp]/t);
          if(s<=t)
            s=z[4*(n0-1)+pp]*(z[4*(n0-2)+2+pp]/(t*(one+sqrt(one+s/t))));
          else
            s=z[4*(n0-1)+pp]*(z[4*(n0-2)+2+pp]/(t+sqrt(t)*sqrt(t+s)));
          t=z[4*(n0-2)+pp]+(s+z[4*(n0-2)+4+pp]);
          z[4*(n0-1)+pp]=z[4*(n0-1)+pp]*(z[4*(n0-2)+pp]/t);
          z[4*(n0-2)+pp]=t;
        }
      }
      z[4*(n0-2)]=z[4*(n0-2)+pp]+sigma;
      z[4*(n0-1)]=z[4*(n0-1)+pp]+sigma;
      n0-=2;
    }

    if (n0==i0)
    {
      // z finished
      // undo shift for last squared diagonal element
      z[4*(n0-1)]=z[4*(n0-1)+pp]+sigma;
      n0--;
      return;
    }
    if (n0<i0)
      // z finished
      return;

    if(pp==2)
      pp=0;

    // reverse the qd-array, if warranted
    if(dmin<=zero || n0<n0in)
    {
      if(cbias*z[4*(i0-1)+pp]<z[4*(n0-1)+pp])
      {
        for(j=0;j<(n0-i0+1)/2;j++)
        {
          t=z[4*(i0+j-1)];
          z[4*(i0+j-1)]=z[4*(n0-j-1)];
          z[4*(n0-j-1)]=t;

          t=z[4*(i0+j-1)+1];
          z[4*(i0+j-1)+1]=z[4*(no-j-1)+1];
          z[4*(n0-j-1)+1]=t;

          t=z[4*(i0+j-1)+2];
          z[4*(i0+j-1)+2]=z[4*(no-j-1)+2];
          z[4*(n0-j-1)+2]=t;

          t=z[4*(i0+j-1)+3];
          z[4*(i0+j-1)+3]=z[4*(no-j-1)+3];
          z[4*(n0-j-1)+3]=t;
        }
        if(n0-i0<=4)
        {
          z[4*(n0-1)+2+pp]=z[4*(i0-1)+2+pp];
          z[4*(n0-1)+3-pp]=z[4*(i0-1)+3-pp];
        }
        dmin2=min(dmin2, z[4*(n0-1)+3-pp]);
        z[4*(n0-1)+2+pp]=min(z[4*(n0-1)+2+pp],
                             min(z[4*(i0-1)+2+pp],z[4*i0+2+pp]));
        z[4*(n0-1)+3-pp]=min(z[4*(n0-1)+3-pp],
                             min(z[4*(i0-1)+3-pp],z[4*i0+3-pp]));
        qmax=max(qmax,max(z[4*(i0-1)+pp],z[4*i0+pp]));
        dmin=-zero;
      }
    }

    // choose a shift
    LASQ4<real_t>(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2,
                  tau, ttype, g);

    // apply dqds until dmin>0
    while (dmin<zero || dmin1<zero)
    {
      LASQ5<real_t>(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2,
                    dn, dn1, dn2, eps);
      ndiv=ndiv+(n0-i0+2);
      iter++;

      // check status
      if(dmin>=zero && dmin1>=zero)
        // convergence
        break;

      if(dmin<zero && dmin1<zero
         && z[4*(n0-2)+2+pp]<tol*(sigma+dn1) && abs(dn)<tol*sigma)
      {
        // convergence hidden by negative dn.
        z[4*(n0-2)-pp+2]=zero;
        dmin=zero;
        break;
      }

      if(dmin<zero)
      {
        // the shift value is too big, select a new tau
        nfail++;
        if(ttype<-22)
        {
          // failed twice
          tau=zero;
        }
        else if(dmin1>zero)
        {
          // late failure, gives excellent shift
          tau=(tau+dmin)*(one-two*eps);
          ttype-=11;
        }
        else
        {
          // early failure, divide previous shift value by 4
          tau=quarter*tau;
          ttype-=12;
        }
      }
      else
      {
        if(dmin!=dmin && tau!=zero)
        {
          // NaN, retry with a zero shift
          tau=zero;
          continue;
        }
        // NaN and zero shift, or possible underflow: call safe dqd
        LASQ6<real_t>(i-, n0, z, pp, dmin, dmin1, dmin2, dn, dn1, dn2);
        ndiv=ndiv+(n0-i0+2);
        iter++;
        tau=zero;
        break;
      }
    }

    // compensated accumulation (gives extra precision)
    if(tau<sigma)
    {
      desig=desig+tau;
      t=sigma+desig;
      desig=desig-(t-sigma);
    }
    else
    {
      t=sigma+tau;
      desig=sigma-(t-tau)+desig;
    }
    sigma=t;
  }
}

#endif
