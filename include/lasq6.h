//
//  lasq6.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq6_h
#define _lasq6_h

/// @file lasq6.h Computes one dqd transform in ping-pong form.

#include <algorithm>
#include "lamch.h"
#include "latl.h"

namespace LATL
{
  /// @brief Computes one dqd (shift equal to zero) transform in ping-pong
  ///  form, with protection against underflow and overflow.
  ///
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] i0 First index.
  /// @param[in] n0 Last index.
  /// @param[in,out] z  Real array of length 4*n. z holds the qd array.
  ///     At output, new emin is stored in z[4*n0-2+pp] to avoid an extra
  ///     argument.
  /// @param[in] pp 0 for ping, 1 for pong.
  /// @param[out] dmin  Minimum value of d.
  /// @param[out] dmin1 Minimum value of d, excluding d[n0-1].
  /// @param[out] dmin2 Minimum value of d, excluding d[n0-1] and d[n0-2].
  /// @param[out] dn    d[n0-1], the last value of d.
  /// @param[out] dnm1  d[n0-2].
  /// @param[out] dnm2  d[n0-3].
  /// @ingroup AUX
  template<typename real_t>
  void LASQ6(const int_t i0, const int_t n0, real_t *z, const int_t pp,
						 real_t &dmin, real_t &dmin1, real_t &dmin2,
						 real_t &dn, real_t &dnm1, real_t &dnm2)
  {
    using std::min;

    if (n0-i0<=1)
      return;

    const real_t zero(0.0);
    const real_t safemin=LAMCH<real_t>('S');

    real_t emin;
    real_t d;
    real_t dmin;
    real_t temp;
    
    emin=z[4*i0+pp]; // why emin=q(i0+1) or qq(i0+1)?
    d=z[4*(i0-1)+pp];
    dmin=d;

    if (pp==0)
    {
      // ping: z[4xj+1] and z[4xj+3] are updated using z[4xj] and z[4xj+2]
      for (j4=4*(i0-1); j4<4*(n0-3); j4+=4)
      {
				z[j4+1]=d+z[j4+2];
				if (z[j4+1]==zero)
				{
					// updated q(j) is zero
					z[j4+3]=zero;
					d=z[j4+4];
					dmin=d;
					emin=zero;
				}
				else if ((safemin*z[j4+4]<z[j4+1]) && (safemin*z[j4+1]<z[j4+4]))
				{
					// no possible under-/overflow: factorize quotient
					temp=z[j4+4]/z[j4+1];
					z[j4+3]=z[j4+2]*temp;
					d=d*temp;
				}
				else
				{
					// avoid under-/overflow
					z[j4+3]=z[j4+4]*(z[j4+2]/z[j4+1]);
					d=z[j4+4]*(d/z[j4+1]);
				}
				dmin=min(dmin,d);
				emin=min(emin,z[j4+3]);
      }
    }
    else
    {
      // pong: z[4xj] and z[4xj+2] are updated using z[4xj+1] and z[4xj+3]
      for (j4=4*(i0-1); j4<4*(n0-3); j4+=4)
      {
				z[j4]=d+z[j4+3];
				if (z[j4]==zero)
				{
					// updated q(j) is zero
					z[j4+2]=zero;
					d=z[j4+5];
					dmin=d;
					emin=zero;
				}
				else if ((safemin*z[j4+5]<z[j4]) && (safemin*z[j4]<z[j4+5]))
				{
					// no possible under-/overflow: factorize quotient
					temp=z[j4+5]/z[j4];
					z[j4+2]=z[j4+3]*temp;
					d=d*temp;
				}
				else
				{
					// avoid under-/overflow
					z[j4+2]=z[j4+5]*(z[j4+3]/z[j4]);
					d=z[j4+5]*(d/z[j4]);
				}
				dmin=min(dmin,d);
				emin=min(emin,z[j4+2]);
      }
    }

    // last two steps are unrolled so as to record last d values
    dnm2=d;
    dmin2=dmin;
    j4=4*(n0-3)+pp; // index of updated values
    j4old=4*(n0-3)-pp+1; // index of old values
    z[j4]=dnm2+z[j4old+2];
    if(z[j4]==zero)
    {
      // updated q(j) is zero
      z[j4+2]=zero;
      dnm1=z[j4old+4];
      dmin=dnm1;
      emin=zero;
    }
    else if((safemin*z[j4old+4]<z[j4]) && (safemin*z[j4]<z[j4old+4]))
    {
      // no possible under-/overflow: factorize quotient
      temp=z[j4old+4]/z[j4];
      z[j4+2]=z[j4old+2]*temp;
      dnm1=dnm2*temp;
    }
    else
    {
      // avoid under-/overflow
      z[j4+2]=z[j4old+4]*(z[j4old+2)/z[j4]);
      dnm1=z[j4old+4]*(dnm2/z[j4]);
    }
    dmin=min(dmin,dmin1);

    dmin1=dmin;
    j4+=4;
    j4old+=4;
    z[j4]=dnm1+z[j4old+2];
    if(z[j4]==zero)
    {
      // updated q(j) is zero
      z[j4+2]=zero;
      dn=z[j4old+4];
      dmin=dn;
      emin=zero;
    }
    else if((safemin*z[j4old+4]<z[j4]) && (safemin*z[j4]<z[j4old+4]))
    {
      // no possible under-/overflow: factorize quotient
      temp=z[j4old+4]/z[j4];
      z[j4+2]=z[j4old+2]*temp;
      dn=dnm1*temp;
    }
    else
    {
      // avoid under-/overflow
      z[j4+2]=z[j4old+4]*(z[j4old+2)/z[j4]);
      dn=z[j4old+4]*(dnm1/z[j4]);
    }
    dmin=min(dmin,dn);

    z[j4+4]=dn;
    z[j4+6]=emin;
  }
}

#endif
