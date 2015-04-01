//
//  lasq5.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq5_h
#define _lasq5_h

/// @file lasq5.h Computes one dqds transform in ping-pong form.

#include <algorithm>
#include "latl.h"

namespace LATL
{
  /// @brief Computes one dqds transform in ping-pong form.
  ///
  /// Assume IEEE-754 arithmetic.
  ///
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] i0 First index.
  /// @param[in] n0 Last index.
  /// @param[in,out] z  Real array of length 4*n. z holds the qd array.
  ///     At output, new emin is stored in z[4*(n0-1)+3-pp] to avoid an extra
	///     argument.
  /// @param[in] pp 0 for ping, 1 for pong.
  /// @param[in,out] tau The shift. It might be set to zero if negligible.
  /// @param[in] sigma The accumulated shift up to this step.
  /// @param[out] dmin  Minimum value of d.
  /// @param[out] dmin1 Minimum value of d, excluding d[n0-1].
  /// @param[out] dmin2 Minimum value of d, excluding d[n0-1] and d[n0-2].
  /// @param[out] dn    d[n0-1], the last value of d.
  /// @param[out] dnm1  d[n0-2].
  /// @param[out] dnm2  d[n0-3].
  /// @param[in] eps    The value of the epsilon used.
  /// @ingroup AUX
  template<typename real_t>
  void LASQ5(const int_t i0, const int_t n0, real_t *z, const int_t pp,
             real_t &tau, const real_t &sigma,
             real_t &dmin, real_t &dmin1, real_t &dmin2,
             real_t &dn, real_t &dnm1, real_t &dnm2, const real_t &eps)
  {
    using std::min;

    if (n0-i0<=1)
      return;

    // z is in qd-format
    // z[4(i-1)]   is qi, z[4(i-1)+1] is qqi,
    // z[4(i-1)+2] is ei, z[4(i-1)+3] is eei

    const real_t zero(0.0);
    const real_t half(0.5);
    const real_t thresh=eps*(sigma+tau);

    real_t emin;
    real_t dmin;
    real_t d;
    real_t temp;

    if(tau<thresh*half)
      // the shift is negligible compared to cumulative shift sigma
      tau=zero;

    emin=z[4*i0+pp]; // emin=q(i0+1) or qq(i0+1)?
    d=z[4*(i0-1)+pp]-tau;
    dmin=d;
    dmin1=-z[4*(i0-1)+pp];

    if(tau!=zero)
    {
      for (int_t j=i0; j<n0-2; j++)
      {
        z[4*(j-1)+1-pp]=d+z[4*(j-1)+2+pp];
        temp=z[4*j+pp]/z[4*(j-1)+1-pp];
        z[4*(j-1)+3-pp]=z[4*(j-1)+2+pp]*temp;
        d=d*temp-tau;
        dmin=min(dmin,d);
        emin=min(emin,z[4*(j-1)+3-pp]);
      }

      // last two steps are unrolled so as to record last d values
      dnm2=d;
      dmin2=dmin;
      z[4*(n0-3)+1-pp]=dnm2+z[4*(n0-3)+2+pp];
      z[4*(n0-3)+3-pp]=z[4*(n0-2)+pp]*(z[4*(n0-3)+2+pp]/z[4*(n0-3)+1-pp]);
      dnm1=z[4*(n0-2)+pp]*(dnm2/z[4*(n0-3)+1-pp])-tau;
      dmin=min(dmin,dmin1);

      dmin1=dmin;
      z[4*(n0-2)+1-pp]=dnm1+z[4*(n0-2)+2+pp];
      z[4*(n0-2)+3-pp]=z[4*(n0-1)+pp]*(z[4*(n0-2)+2+pp]/z[4*(n0-2)+1-pp]);
      dn=z[4*(n0-1)+pp]*(dnm2/z[4*(n0-2)+1-pp])-tau;
      dmin=min(dmin,dn);
    }
    else
    {
      // set d's to zero when negligible compared to sigma
      for (int_t j=i0; j<n0-2; j++)
      {
        z[4*(j-1)+1-pp]=d+z[4*(j-1)+2+pp];
        temp=z[4*j+pp]/z[4*(j-1)+1-pp];
        z[4*(j-1)+3-pp]=z[4*(j-1)+2+pp]*temp;
        d=d*temp;
        if(d<thresh)
          d=zero;
        dmin=min(dmin,d);
        emin=min(emin,z[4*(j-1)+3-pp]);
      }

      // last two steps are unrolled so as to record last d values
      dnm2=d;
      dmin2=dmin;
      z[4*(n0-3)+1-pp]=dnm2+z[4*(n0-3)+2+pp];
      if(dnm2<zero)
        return;
      z[4*(n0-3)+3-pp]=z[4*(n0-2)+pp]*(z[4*(n0-3)+2+pp]/z[4*(n0-3)+1-pp]);
      dnm1=z[4*(n0-2)+pp]*(dnm2/z[4*(n0-3)+1-pp]);
      dmin=min(dmin,dmin1);

      dmin1=dmin;
      z[4*(n0-2)+1-pp]=dnm1+z[4*(n0-2)+2+pp];
      if(dnm1<zero)
        return;
      z[4*(n0-2)+3-pp]=z[4*(n0-1)+pp]*(z[4*(n0-2)+2+pp]/z[4*(n0-2)+1-pp]);
      dn=z[4*(n0-1)+pp]*(dnm2/z[4*(n0-2)+1-pp]);
      dmin=min(dmin,dn);
    }

    z[4*(n0-1)+pp]=dn;
    z[4*(n0-1)+2+pp]=emin;
  }
}

#endif
