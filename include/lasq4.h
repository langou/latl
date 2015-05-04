//
//  lasq4.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq4_h
#define _lasq4_h

/// @file lasq4.h Computes an approximation to the smallest eigenvalue using
///   values of d from the previous transform.

#include <algorithm>
#include <cmath>
#include "latl.h"

namespace LATL
{
  /// @brief Computes an approximation tau to the smallest eigenvalue using
  ///   values of d from the previous transform.
  ///
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] i0   First index.
  /// @param[in] n0   Last index.
  /// @param[in] z    Real array of length 4*n. z holds the qd array. emin is
  ///   stored in z[4*n0-1] to avoid an extra argument.
  /// @param[in] pp 0 for ping, 1 for pong.
  /// @param[in] n0in  The value of n0 at start of eigtest.
  /// @param[in] dmin  Minimum value of d.
  /// @param[in] dmin1 Minimum value of d, excluding d[n0-1].
  /// @param[in] dmin2 Minimum value of d, excluding d[n0-1] and d[n0-2].
  /// @param[in] dn    d[n0-1], the last value of d.
  /// @param[in] dn1   d[n0-2].
  /// @param[in] dn2   d[n0-3].
  /// @param[out] tau  The shift.
  /// @param[in,out] ttype Shift type.
  ///   On entry, tttype<-11 signals a failure.
  ///   Especially, ttype=-18 means first failure in case 6
  ///   On exit,
  ///   ttype=-1, new qd segment
  ///   ttype=-2, asymptotic case (shift estimated with the Rayleigh quotient)
  ///   ttype=-3, asymptotic case (shift estimated with Gershgorin's circles)
  ///   ttype=-4, not quite asymptotic case
  ///   ttype=-5, dmin is d[n0-2]
  ///   ttype=-6, typical situation in early stages
  ///   ttype=-7, one eigenvalue deflated (Rayleigh quotient with gap)
  ///   ttype=-8, one eigenvalue deflated (Rayleigh quotient)
  ///   ttype=-9, one eigenvalue deflated and non-asymptotic case
  ///   ttype=-10, two eigenvalues deflated (Rayleigh quotient)
  ///   ttype=-11, two eigenvalues deflated
  ///   ttype=-12, more than two eigenvalues deflated
  /// @param[in,out] g     On entry, last fraction of dmin used in case 6.
  ///   On exit, new fraction of dmin.
  /// @ingroup AUX
  template<typename real_t>
  void LASQ4(int_t i0, int_t n0, real_t *z, int_t pp, int_t n0in,
             real_t &dmin, real_t &dmin1, real_t &dmin2,
             real_t &dn, real_t &dn1, real_t &dn2,
             real_t &tau, int_t &ttype, real_t &g)
  {
    using std::min;
    using std::max;
    using std::sqrt;
    const real_t zero(0.0);

    if(dmin<=zero)
    {
    // case 1, new qd segment
      // a negative dmin forces the shift to take its absolute value
      tau=-dmin;
      ttype=-1;
      return;
    }

    const real_t one(1.0);
    const real_t two(2.0);
    const real_t three(3.0);
    const real_t hundred(100.0);
    const real_t half(0.5);
    const real_t third(one/three);
    const real_t quarter(0.25);
    const real_t cnst2(1.01);
    const real_t cnst3(1.05);

    // z is in qd-format
    // z[4(i-1)]   is qi, z[4(i-1)+1] is qqi,
    // z[4(i-1)+2] is ei, z[4(i-1)+3] is eei

    if(n0in==n0)
    {
      // no eigenvalue deflated
      if(dmin==dn || dmin==dn1)
      {
        if(dmin==dn && dmin1==dn1)
        {
          // asymptotic cases
          real_t a2, b1, b2;
          real_t gap1, gap2;

          b1=sqrt(z[4*(n0-1)+pp])*sqrt(z[4*(n0-2)+2+pp]);
          b2=sqrt(z[4*(n0-2)+pp])*sqrt(z[4*(n0-3)+2+pp]);
          a2=z[4*(n0-2)+pp]+z[4*(n0-2)+2+pp];
          gap2=dmin2-a2-dmin2*quarter;
          if(gap2>zero && gap2>b2)
            gap1=a2-dn-(b2/gap2)*b2;
          else
            gap1=a2-dn-(b1+b2);
          if(gap1>zero && gap1>b1)
          {
            // case 2, estimate lower bound with Rayleigh quotient
            tau=max(dn-(b1/gap1)*b1, half*dmin);
            ttype=-2;
            return;
          }
          else
          {
            // case 3, estimate lower bound with Gershgorin circles
            // warning: LAPACK 3.4.2 code does not match the algorithm
            //    described in the article referenced above.
            //    Here, we follow the LAPACK code and set tau to
            //    max(dn/3, dn-b1) when a2<=b1+b2.
            tau=zero;
            if(dn>b1)
              tau=dn-b1;
            if(a2>(b1+b2))
              tau=min(tau,a2-(b1+b2));
            tau=max(tau, third*dmin);
            ttype=-3;
            return;
          }
        }
        else
        {
          // case 4, not quite asymptotic
          tau=quarter*dmin; //FIXME: LAPACK uses Rayleigh quotient
          ttype=-4;
          return;
        }
      }
      else if(dmin==dn2)
      {
        // case 5
        tau=quarter*dmin; // FIXME:LAPACK uses twisted factorization
        ttype=-5;
        return;
      }
      else
      {
        // case 6, no information to guide us
        // warning: LAPACK 3.4.2 code does not match the algorithm
        //   described in the article referenced above.
        //   Here, we follow the LAPACK code and set tau to
        //   (1/3+2*f/3)*tau for the last step.
        if(ttype==-6)
          // last step
          g=g+third*(one-g);
        else if(ttype==-18)
          // first failure
          g=quarter*third;
        else
          g=quarter;
        tau=g*dmin;
        ttype=-6;
        return;
      }
    }
    else if(n0in==n0+1)
    {
      // one eigenvalue just deflated
      if(dmin1==dn1 && dmin2==dn2)
      {
        // asymptotic cases
        tau=third*dmin1;
        if(z[4*(n0-2)+2+pp]>z[4*(n0-2)+pp])
        {
          ttype=-7;
          return;
        }
        real_t zi, sum, tmp;
        zi=z[4*(n0-2)+2+pp]/z[4*(n0-2)+pp];
        sum=zi;
        if(sum!=zero)
        {
          for(int_t j=n0-2; j>=i0; j--)
          {
            tmp=zi;
            if(z[4*(j-1)+2+pp]>z[4*(j-1)+pp])
            {
              ttype=-7;
              return;
            }
            zi=zi*z[4*(j-1)+2+pp]/z[4*(j-1)+pp];
            sum=sum+zi;
            if(hundred*max(zi,tmp)<sum)
              break;
          }
        }
        real_t gap;
        sum=cnst3*sum;
        tmp=dmin1/(one+sum);
        gap=half*dmin2-tmp;
        sum=sqrt(sum);
        if(gap>zero && gap>sum*tmp)
        {
          // case 7, use Rayleigh quotient and gap
          tau=max(tau, tmp*(one-cnst2*tmp*(sum/gap)*sum));
          ttype=-7;
          return;
        }
        else
        {
          // case 8, use Rayleigh quotient
          tau=max(tau, tmp*(one-cnst2*sum));
          ttype=-8;
          return;
        }
      }
      else
      {
        // case 9, non-asymptotic case
        tau=quarter*dmin1;
        if(dmin1==dn1)
          tau=half*dmin1;
        ttype=-9;
        return;
      }
    }
    else if(n0in==n0+2)
    {
      // two eigenvalues deflated
      if(dmin2==dn2 && two*z[4*(n0-2)+2+pp]<z[4*(n0-2)+pp])
      {
        // case 10, use Rayleigh quotient
        tau=third*dmin2;
        ttype=-10;
        if(z[4*(n0-2)+2+pp]>z[4*(n0-2)+pp])
          return;

        real_t zi, sum, tmp;
        zi=z[4*(n0-2)+2+pp]/z[4*(n0-1)+pp];
        sum=zi;
        if(sum!=zero)
        {
          for(int_t j=n0; j>i0; j--)
          {
            if(z[4*(j-1)+2+pp]>z[4*(j-1)+pp])
              return;
            zi=zi*z[4*(j-1)+2+pp]/z[4*(j-1)+pp];
            sum=sum+zi;
            if(hundred*zi<sum)
              break;
          }
        }
        real_t gap;
        sum=cnst3*sum;
        tmp=dmin2/(one+sum);
        sum=sqrt(sum);
        gap=z[4*(n0-2)+pp]+z[4*(n0-3)+2+pp]
          -sqrt(z[4*(n0-3)+pp])*sqrt(z[4*(n0-3)+2+pp])-tmp;
        if(gap>zero && gap>sum*tmp)
          tau=max(tau, tmp*(one-cnst2*tmp*(sum/gap)*sum));
        else
          tau=max(tau, tmp*(one-cnst2*sum));
        return;
      }
      else
      {
        // case 11, non-asymptotic case
        tau=quarter*dmin2;
        ttype=-11;
        return;
      }
    }
    else
    {
      // case 12, more than two eigenvalues deflated
      tau=zero;
      ttype=-12;
      return;
    }
  }
}

#endif
