//
//  lasq2.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 3/17/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasq2_h
#define _lasq2_h

/// @file lasq2.h Computes all eigenvalues of the symmetric positive definite
///  tridiagonal matrix associated with the qd array z to high relative
///  accuracy.

#include <algorithm>
#include <cmath>
#include "lasq3.h"
#include "lasrt.h"
#include "lamch.h"
#include "latl.h"

namespace LATL
{
  /// @brief Computes all the eigenvalues of the symmetric positive definite
  ///  tridiagonal matrix associated with the qd array z to high relative
  ///  accuracy, in the absence of denormalization, underflow, and overflow.
  ///
  /// To see the relation of z to the tridiagonal matrix, let L be a unit
  /// lower bidiagonal with subdiagonals z[2,4,6,...] and let U be an upper
  /// bidiagonal matrix with 1's above and diagonal z[1,3,5,...].
  /// The related tridiagonal is the symmetric tridiagonal similar to L*U.
  ///
  /// See "An implementation of the dqds algorithm (positive case)" by B. N.
  /// Parlett and O. A. Marques, LAPACK Working Note 155.
  ///
  /// @tparam real_t Floating point type.
  /// @param[in] n      The number of rows and columns in the matrix. n>=0.
  /// @param[in,out] z  Array of length 4*n.
  ///   On entry, z holds the qd array.
  ///   On successful exit, entries 1 to n hold the eigenvalues in decreasing
  ///   order, z[2*n] holds the trace, and z[2*n+1] holds the sum of the
  ///   eigenvalues (IS IT THE TRACE,ISN'T IT?). If n>2, then z[2*n+2] holds
  ///   the iteration count, z[2*n+3] holds ndivs/nin^2, and z[2*n+4] holds
  ///   the percentage of shifts that failed.
  ///   On exit with error,
  /// @return 0 if success.
  /// @return -i if the ith argument is a scalar and had an illegal value.
  /// @return -(100*i+j) if the i-th argument is an array and the j-entry
  ///    (0-based) had an illegal value.
  /// @return 1 if a split was marked by a positive value in e
  /// @return 2 if current block of z is not diagonalized after 100*n
  ///   iterations in inner while loop. On exit, z holds a qd array with the
  ///   same eigenvalues as the given z.
  /// @return 3 if the termination criterion of outer while loop is not met
  ///   (the program created more than n unreduced blocks).
  /// @ingroup AUX

  template<typename real_t>
  int_t LASQ2(int_t n, real_t *z)
  {
    const real_t zero(0.0);

    using std::abs;
    using std::sqrt;
    using std::max;
    using std::min;

    if(n<0)
      return -1;

    if(n==0)
      return;
    if(n==1)
    {
      // 1-by-1 case
      if(z[0]<zero)
        return -201;
      return 0;
    }

    const real_t one(1.0);
    const real_t hundred(100.0);
    const real_t half(0.5);
    const real_t eps=LAMCH<real_t>('P');
    const real_t safemin=LAMCH<real_t>('S');
    const real_t tol=eps*hundred;
    const real_t tol2=tol*tol;

    if(n==2)
    {
      // 2-by-2 case
      if(z[1]<zero || z[2]<zero)
        return -2;

      real_t d;

      if(z[2]>z[0])
      {
        d=z[2];
        z[2]=z[0];
        z[0]=d;
      }
      z[4]=z[0]+z[1]+z[2]; // trace of LU before transformation
      if(z[1]>z[2]*tol2)
      {
        // non-negligible e1
        // compute two eigenvalues following Rutishauser's algorithm
        d=half*((z[0]-z[2])+z[1]);
        if(d!=zero)
        {
          real_t s;
          s=z[2]*(z[1]/d);
          if(s<=t)
            s=z[2]*(z[1]/(t*(one+sqrt(one+s/d))));
          else
            s=z[2]*(z[1]/(d+sqrt(d)*sqrt(d+s)));
          t=z[0]+(s+z[1]);
          z[2]=z[2]*(z[0]/d);
          z[0]=d;
        }
      }
      z[1]=z[2];
      z[5]=z[0]+z[1]; // trace after diagonalisation
      return 0;
    }

    // check for negative data and compute sums of q's and e's
    real_t emin;
    real_t qmax;
    real_t zmax;
    real_t d;
    real_t e;

    z[2*n-1]=zero;
    emin=z[1];
    qmax=zero;
    zmax=zero;
    d=zero;
    e=zero;

    for(int_t k=0; k<2*(n-1); k+=2)
    {
      if(z[k]<zero)
        return -(200+k+1);
      if(z[k+1]<zero)
        return -(200+k+2);
      d=d+z[k];
      e=e+z[k+1];
      qmax=max(qmax, z[k]);
      emin=min(emin,z[k+1]);
      zmax=max(zmax,max(z[k],z[k+1]));
    }
    if(z[2*n-2]<zero)
      return -(200+2*n-1);
    d=d+z[2*n-2];
    qmax=max(qmax, z[2*n-1]);
    zmax=max(zmax, z[2*n-1]);

    // check for diagonality
    if(e==zero)
    {
      for(int_t k=1; k<n; k++)
        z[k]=z[2*k];
      LASRT<real_t>('D',n,z);
      z[2*n-2]=d;
      return 0;
    }

    // check for zero data
    real_t trace=d+e;
    if(trace==zero)
    {
      z[2*n-2]=zero;
      return 0;
    }

    // rearrange for locality in a qd array z
    // z[4(i-1)] is qi, z[4(i-1)+1] is qqi,
    // z[4(i-1)+2] is ei, and z[4(i-1)+3] is eei
    for(k=2*n-1, k>=0; k-=2)
    {
      z[2*k+3]=zero;
      z[2*k+2]=z[k+1];
      z[2*k+1]=zero;
      z[2*k]=z[k];
    }

    int_t i0=1;
    int_t n0=n;

    const real_t cbias(1.5);
    if(cbias*z[4*(i0-1)]<z[4*(n0-1)])
    {
      // reverse the qd-array
      for(int_t k=0; k<n0-i0; k++)
      {
        d=z[4*(no-k-1)];
        z[4*(n0-k-1)]=z[4*(i0+k-1)];
        z[4*(i0+k-1)]=d;
        d=z[4*(no-k-1)+2];
        z[4*(n0-k-1)+2]=z[4*(i0+k-1)+2];
        z[4*(i0+k-1)+2]=d;
      }
    }

    // initial split via dqd and Li's test
    // in a ping pass (pp=0), (qqi,eei) are computed from (qi,ei)
    // in a pong pass (pp=1), (qi,ei) are computed from (qqi,eei)
    real_t tmp;
    int_t pp=0;
    for(int_t k=1; k<3; k++)
    {
      d=z[4*(n0-1)+pp];
      for(int_t i=n0-1; i>=i0-1; i--)
      {
        if(z[4*i+2+pp]<=tol2*d)
        {
          z[4*i+2+pp]=-zero;
          d=z[4*i+pp];
        }
        else
          d=z[4*i+pp]*(d/(d+z[4*i+2+pp]));
      }

      // safe dqd
      emin=z[4*(i0-1)+2+pp];
      d=z[4*(i0-1)+pp];
      for(int_t i=i0-1; i<n0; i++)
      {
        z[4*i+1-pp]=d+z[4*i+2+pp];
        if(z[4*i+2+pp]<=tol2*d)
        {
          z[4*i+2+pp]=-zero;
          z[4*i+pp]=d;
          z[4*i+3-pp]=zero;
          d=z[4*i+4+pp];
        }
        else if(safemin*z[4*i+4+pp]<z[4*i+1-pp]
                && safemin*z[4*i+1-pp]<z[4*i+4-pp])
        {
          tmp=z[4*i+4+pp]/z[4*i+1-pp];
          z[4*i+2+pp]=z[4*i+3-pp]*tmp;
          d=d*tmp;
        }
        else
        {
          z[4*i+2+pp]=z[4*i+4+pp]*(z[4*i+3-pp]/z[4*i+1-pp]);
          d=z[4*i+4+pp]*(d/z[4*i+1-pp]);
        }
        emin=min(emin,z[4*i+2+pp]);
      }
      z[4*(n0-1)+1-pp]=d;

      qmax=z[4*(i0-1)+pp];
      for(int_t i=i0; i<n0; i++)
        qmax=max(qmax,z[4*i+pp]);

      pp=1-pp;
    }

    const real_t four(4.0);
    int_t ttype=0;
    real_t dmin1=zero;
    real_t dmin2=zero;
    real_t dn=zero;
    real_t dn1=zero;
    real_t dn2=zero;
    real_t g=zero;
    real_t tau=zero;
    real_t sigma;
    real_t desig;

    int_t j;
    int_t iter=2;
    int_t nfail=0;
    int_t ndiv=2*(n0-i0);
    real_t emax;

    for(int_t i=1; i<n+2; i++)
    {
      if(n0<1)
      {
        // move q's to the front
        for(j=1;j<n;j++)
          z[j]=z[4*j];
        // sort and compute eigenvalues
        LASRT<real_t>('D',n,z);

        // store trace, sum(eigenvalues) and information on performance
        z[2*n]=trace;
        z[2*n+1]=zero;
        for(j=n-1; j>=0; j++)
          z[2*n+1]=z[2*n+1]+z[j];
        z[2*n+2]=real_t(iter);
        z[2*n+3]=real_t(ndiv)/real_t(n*n);
        z[2*n+4]=hundred*real_t(nfail)/real_t(iter);
        return 0;
      }

      desig=zero;
      if(n0=n)
        sigma=zero;
      else
        // e(n0) holds the value of sigma when submatrix in i0:n0
        // splits from the rest of the array, but is negated
        sigma=-z[4*(n0-1)+2];
      if(sigma<zero)
        return 1;

      // find qmax and emin
      // find Gershgorin-type bound if q's are much greater than e's
      emax=zero;
      if(n0>i0)
        emin=abs(z[4*(n0-2)+2]);
      else
        emin=zero;
      qmin=z[4*(n0-1)];
      qmax=qmin;
      for(j=n0-2; j>0; j--)
      {
        // find top index i0 of the last unreduced submatrix
        if(z[4*j+2]<=zero)
          break;

        if(qmin>=four*emax)
        {
          qmin=min(qmin,z[4*j+4]);
          emax=max(emax,z[4*j+2]);
        }
        qmax=max(qmax,z[4*j]+z[4*j+2]);
        emin=min(emin,z[4*j+2]);
      }
      i0=j+2;
      pp=0;

      if(n0-i0>1)
      {
        int_t kmin;
        real_t dee;
        real_t deemin;

        dee=z[4*(i0-1)];
        deemin=dee;
        kmin=i0;
        for(j=i0; j<n0-3; j++)
        {
          dee=z[4*(j-1)]*(dee/(dee+z[4*j+2]));
          if(dee<=deemin)
          {
            deemin=dee;
            kmin=j+1;
          }
        }
        if((kmin-i0)*2<n0-kmin && deemin<=half*z[4*(n0-1)])
        {
          // flip segment i0:n0
          pp=2;
          for(j=0; j<(n0-i0+1)/2; j++)
          {
            tmp=z[4*(no-j-1)];
            z[4*(n0-j-1)]=z[4*(i0+j-1)];
            z[4*(i0+j-1)]=tmp;
            tmp=z[4*(no-j-1)+1];
            z[4*(n0-j-1)+1]=z[4*(i0+j-1)+1];
            z[4*(i0+j-1)+1]=tmp;
            tmp=z[4*(no-j-1)+2];
            z[4*(n0-j-1)+2]=z[4*(i0+j-1)+2];
            z[4*(i0+j-1)+2]=tmp;
            tmp=z[4*(no-j-1)+3];
            z[4*(n0-j-1)+3]=z[4*(i0+j-1)+3];
            z[4*(i0+j-1)+3]=tmp;
          }
        }
      }

      // put -(initial shift) to dmin
      dmin=-max(zero,qmin-two*sqrt(qmin)*sqrt(emax));

      // now i0:n0 is unreduced
      int_t nbig=100*(n0-i0+1);
      for(int_t it=0; it<nbig; it++)
      {
        if(i0>n0)
          break;

        // apply dqds to z
        LASQ3<real_t>(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv,ttype,
                      dmin1,dmin2,dn,dn1,dn2,g,tau);
        pp=1-pp;

        if(pp==0 && n0-i0>3)
        {
          if(z[4*(n0-1)+3]<=tol2*qmax
             || z[4*(n0-1)+2]<=tol2*sigma)
          {
            // emin is very small, check for splits
            int_t split=i0-1;
            qmax=z[4*(i0-1)];
            emin=z[4*(i0-1)+2];
            tmp=z[4*(i0-1)+3]; //old emin
            for(j=i0; j<n0-3; j++)
            {
              if(z[4*(j-1)+3]<=tol2*z[4*(j-1)]
                 || z[4*(j-1)+2]<=tol2*sigma)
              {
                z[4*(j-1)+2]=-sigma;
                split=j;
                qmax=zero;
                emin=z[4*j+2];
                tmp=z[4*j+3];
              }
              else
              {
                qmax=max(qmax,z[4*j]);
                emin=min(emin,z[4*(j-1)+2]);
                tmp=min(tmp,z[4*(j-1)+3]);
              }
            }
            z[4*(n0-1)+2]=emin;
            z[4*(n0-1)+3]=tmp;
            i0=split+1;
          }
        }

        if(i0<=n0)
        {
          // the number of iterations have exceeded the maximum
          int_t i1=i0;
          int_t n1=n0;
          real_t tmpe;

          // restore the shift sigma
          while(i1>1)
          {
            tmp=z[4*(i0-1)];
            z[4*(i0-1)]=z[4*(i0-1)]+sigma;
            for(j=i0-1;j<n0;j++)
            {
              tmpe=z[4*j+2];
              z[4*j+2]=z[4*j+2]*(tmp/z[4*j]);
              tmp=z[4*j+4];
              z[4*j+4]=z[4*j+4]+sigma+tmpe-z[4*j+2];
            }
            if(i1>1)
            {
              n1=i1-1;
              while(i1>=2 && z[4*(i1-2)+2])
                i1--;
              sigma=-z[4*(n1-1)+2];
            }
          }

          // place the new d's and e's in a qd-array
          for(j=0;j<n;j++)
          {
            z[2*j]=z[4*j];
            if(j<n0)
              // copy e's from the unfinished block
              z[2*j+1]=z[4*j+2];
            else
              // other data may have been stored in e's places
              z[2*j+1]=zero;
          }
          return 2;
        }
      }
    }
    return 3;
  }
}

#endif
