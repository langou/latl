//
//  lasr.h
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 1/7/15.
//  Copyright (c) 2015 University of Colorado Denver. All right reserved.
//

#ifndef _lasr_h
#define _lasr_h

/// @file lasr.h Applies a sequence of plane rotations to a general
///   rectangular matrix.

#include <cctype>
#include "latl.h"

namespace LATL
{
  /// @brief The real version applies a sequence of plane rotations to a real
  ///   matrix A, from either the left or the right.
  ///
  /// When side='L', the transformation takes the form
  ///
  ///   A = P * A
  ///
  /// and when side='R', the transformation takes the form
  ///
  ///   A = A * P'
  ///
  /// where P is an orthogonal matrix consisting of a sequence of z plane
  /// rotations, with z=m when side='L' and z=n when side='R'.
  ///
  /// When direct='F' (Forward sequence), then
  ///
  ///   P = P(z-1) * ... * P(2) * P(1)
  ///
  /// and when direct='B' (Backward sequence), then
  ///
  ///   P = P(1) * P(2) * ... * P(z-1)
  ///
  /// where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
  ///
  ///   R(k) = ( c(k)  s(k))
  ///          (-s(k)  c(k)).
  ///
  /// When pivot='V' (Variable pivot), the rotation is performed for the
  /// plane (k, k+1), i.e. P(k) has the form
  ///
  ///   P(k) = ( 1                             )
  ///          (   ...                         )
  ///          (       1                       )
  ///          (           c(k)  s(k)          )
  ///          (          -s(k)  c(k)          )
  ///          (                       1       )
  ///          (                         ...   )
  ///          (                             1 )
  ///
  /// where R(k) appears as a rank-2 modification to the identity matrix in
  /// rows and columns k and k+1.
  ///
  /// When pivot='T' (Top pivot), the rotation is performed for the
  /// plane (1, k+1), i.e. P(k) has the form
  ///
  ///   P(k) = (  c(k)           s(k)          )
  ///          (       1                       )
  ///          (         ...                   )
  ///          (             1                 )
  ///          ( -s(k)           c(k)          )
  ///          (                       1       )
  ///          (                         ...   )
  ///          (                             1 )
  ///
  /// where R(k) appears in rows and columns 1 and k+1.
  ///
  /// When pivot='B' (Bottom pivot), the rotation is performed for the
  /// plane (k, z), i.e. P(k) has the form
  ///
  ///   P(k) = ( 1                             )
  ///          (   ...                         )
  ///          (       1                       )
  ///          (           c(k)           s(k) )
  ///          (                  1            )
  ///          (                    ...        )
  ///          (                        1      )
  ///          (          -s(k)           c(k) )
  ///
  /// where R(k) appears in rows and columns k and z.
  /// The rotations are performed without forming P(k) explicitely.
  ///
  /// @return 0 if success.
  /// @return -i if the ith argument is invalid.
  /// @tparam real_t Floating point type.
  /// @param[in] side   Specifies whether the plane rotation matrix P is
  ///   applied to A on the left or on the right.
  ///   If side='L' or 'l' (Left), compute A = P * A.
  ///   If side='R' or 'r' (Right), compute A = A * P'.
  /// @param[in] pivot  Specifies the plane for which P(k) is a plane
  ///   rotation matrix.
  ///   If pivot='V' or 'v' (Variable), the plane (k, k+1).
  ///   If pivot='T' or 't' (Top), the plane (1, k+1).
  ///   If pivot='B' or 'b' (Bottom), the plane (k, z)
  /// @param[in] direct Specifes whether P is a forward or backward sequence
  ///   of plane rotations.
  ///   If direct='F' or 'f' (Forward), P = P(z-1) * ... * P(2) * P(1).
  ///   If direct='B' or 'b' (Backward), P = P(1) * P(2) * ... * P(z-1).
  /// @param[in] m      The number of rows of the matrix A. m>=0.
  /// @param[in] n      The number of columns of the matrix A. n>=0.
  /// @param[in] c      The cosines c(k) of the plane rotations.
  ///   Array of dimension m-1 if side='L' or n-1 if side='R'.
  /// @param[in] s      The sines s(k) of the plane rotation.
  ///   Array of dimension m-1 if side='L' or n-1 if side='R'.
  /// @param[in,out] A  Pointer to a real m-by-n matrix A.
  ///   On exit, A is overwritten by P * A if side='R' or A * P' if side='L'.
  /// @param[in] ldA    The leading dimension of the array A. ldA>=max(1, m).
  /// @ingroup AUX

  template<typename real_t>
  int_t LASR(char side, char pivot, char direct, int_t m, int_t n,
             real_t *c, real_t *s, real_t *A, int_t ldA)
  {
    using std::toupper;
    const real_t one(1.0);
    const real_t zero(0.0);
    int_t i,j;
    real_t temp, ctemp, stemp;

    side=toupper(side);
    pivot=toupper(pivot);
    direct=touper(direct);

    if((side!='L') && (side!='R'))
      return -1;
    else if((pivot!='V') && (pivot!='T') && (pivot!='B'))
      return -2;
    else if((direct!='F') && (direct!='B'))
      return -3;
    else if(m<0)
      return -4;
    else if(n<0)
      return -5;
    else if(ldA<m || ldA<1)
      return -9;

    if(m==0 || n==0)
      // quick return
      return 0;

    if(side=='L')
    {
      // Form P * A
      if(pivot=='V')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (k, k+1)
          for(j=0;j<m-1;j++)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+1+i*ldA];
                A[j+1+i*ldA]=ctemp*temp-stemp*A[j+i*ldA];
                A[j+i*ldA]=stemp*temp+ctemp*A[j+i*ldA];
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (k, k+1)
          for(j=m-2;j>=0;j--)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+1+i*ldA];
                A[j+1+i*ldA]=ctemp*temp-stemp*A[j+i*ldA];
                A[j+i*ldA]=stemp*temp+ctemp*A[j+i*ldA];
              }
            }
          }
        }
      }
      else if(pivot=='T')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (1, k+1)
          for(j=1;j<m;j++)
          {
            ctemp=c[j-1];
            stemp=s[j-1];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+i*ldA];
                A[j+i*ldA]=ctemp*temp-stemp*A[i*ldA];
                A[i*ldA]=stemp*temp+ctemp*A[i*ldA];
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (1, k+1)
          for(j=m-1;j>0;j--)
          {
            ctemp=c[j-1];
            stemp=s[j-1];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+i*ldA];
                A[j+i*ldA]=ctemp*temp-stemp*A[i*ldA];
                A[i*ldA]=stemp*temp+ctemp*A[i*ldA];
              }
            }
          }
        }
      }
      else if(pivot=='B')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (k, m)
          for(j=0;j<m-1;j++)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+i*ldA];
                A[j+i*ldA]=stemp*A[m-1+i*ldA]+ctemp*temp;
                A[m-1+i*ldA]=ctemp*A[m-1+i*ldA]-stemp*temp;
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (k, m)
          for(j=m-2;j>=0;j--)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<n;i++)
              {
                temp=A[j+i*ldA];
                A[j+i*ldA]=stemp*A[m-1+i*ldA]+ctemp*temp;
                A[m-1+i*ldA]=ctemp*A[m-1+i*ldA]-stemp*temp;
              }
            }
          }
        }
      }
    }
    else if(side=='R')
    {
      // Form A * P'
      if(pivot=='V')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (k, k+1)
          for(j=0;j<n-1;j++)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+(j+1)*ldA];
                A[i+(j+1)*ldA]=ctemp*temp-stemp*A[i+j*ldA];
                A[i+j*ldA]=stemp*temp+ctemp*A[i+j*ldA];
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (k, k+1)
          for(j=n-2;j>=0;j--)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+(j+1)*ldA];
                A[i+(j+1)*ldA]=ctemp*temp-stemp*A[i+j*ldA];
                A[i+j*ldA]=stemp*temp+ctemp*A[i+j*ldA];
              }
            }
          }
        }
      }
      else if(pivot=='T')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (1, k+1)
          for(j=1;j<n;j++)
          {
            ctemp=c[j-1];
            stemp=s[j-1];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+j*ldA];
                A[i+j*ldA]=ctemp*temp-stemp*A[i];
                A[i]=stemp*temp+ctemp*A[i];
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (1, k+1)
          for(j=n-1;j>0;j--)
          {
            ctemp=c[j-1];
            stemp=s[j-1];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+j*ldA];
                A[i+j*ldA]=ctemp*temp-stemp*A[i];
                A[i]=stemp*temp+ctemp*A[i];
              }
            }
          }
        }
      }
      else if(pivot=='B')
      {
        if(direct=='F')
        {
          // Forward sequence in plane (k, n)
          for(j=0;j<n-1;j++)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+j*ldA];
                A[i+j*ldA]=stemp*A[i+(n-1)*ldA]+ctemp*temp;
                A[i+(n-1)*ldA]=ctemp*A[i+(n-1)*ldA]-stemp*temp;
              }
            }
          }
        }
        else if(direct=='B')
        {
          // Backward sequence in plane (k, n)
          for(j=n-2;j>=0;j--)
          {
            ctemp=c[j];
            stemp=s[j];
            if((ctemp!=one)||(stemp!=zero))
            {
              for(i=0;i<m;i++)
              {
                temp=A[i+j*ldA];
                A[i+j*ldA]=stemp*A[i+(n-1)*ldA]+ctemp*temp;
                A[i+(n-1)*ldA]=ctemp*A[i+(n-1)*ldA]-stemp*temp;
              }
            }
          }
        }
      }
    }
  }
}

#endif
