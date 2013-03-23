//
//  rotmg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _rotmg_h
#define _rotmg_h

/// @file rotmg.h Constructs modified Givens transformation matrix.

#include <cmath>

namespace LATL
{
   /// @brief Constructs the modified Givens transformation matrix H.
   ///
   /// The Givens rotation matrix H is constructed that zeros the second
   /// component of the vector
   ///
   ///        (sqrt(d1)*x1,sqrt(d2)*y1).
   /// The value of flag indicats the form of H as follows:
   ///
   ///          flag=-1          flag=0          flag=1         flag=-2
   ///        ------------    ------------    -------------    ---------
   ///        [ h11  h12 ]    [  1   h12 ]    [ h11    1  ]    [ 1   0 ]
   ///        [ h21  h22 ]    [ h21   1  ]    [ -1    h22 ]    [ 0   1 ]
   /// and is used in the application routine LATL::ROTM.
   /// @tparam real_t Floating point type.
   /// @param d1 Real scalar.
   /// @param d2 Real scalar.
   /// @param x1 Real scalar.
   /// @param y1 Real scalar.
   /// @param param Pointer to flag and 2-by-2 Givens transformation matrix.
   /// @ingroup ROT
   
   template <typename real_t>
   int ROTMG(real_t &d1, real_t &d2, real_t &x1, real_t y1, real_t *param)
   {
      using std::abs;

      const real_t zero(0.0);
      const real_t one(1.0);
      const real_t two(2.0);
      const real_t gam(4096.0);
      const real_t gamsq(16777216.0);
      const real_t rgamsq(5.9604645e-08);
      real_t u=zero;
      real_t p1=zero;
      real_t p2=zero;
      real_t q1=zero;
      real_t q2=zero;
      real_t h11=zero;
      real_t h12=zero;
      real_t h21=zero;
      real_t h22=zero;
      real_t temp=zero;
      real_t flag=zero;
      
      if(d1<zero)
      {
         flag = -one;
         h11 = zero;
         h12 = zero;
         h21 = zero;
         h22 = zero;
         d1 = zero;
         d2 = zero;
         x1 = zero;
      }
      else
      {
         p2 = d2*y1;
         if(p2==zero)
         {
            flag=-two;
            param[0]=flag;
            return 0;
         }
         p1 = d1*x1;
         q2 = p2*y1;
         q1 = p1*x1;
         if(abs(q1)>abs(q2))
         {
            h21 = -y1/x1;
            h12 = p2/p1;
            u = one - h12*h21;
            if(u > zero)
            {
               flag = zero;
               d1 = d1/u;
               d2 = d2/u;
               x1 = x1*u;
            }
         }
         else
         {
            if(q2 < zero)
            {
               flag = -one;
               h11 = zero;
               h12 = zero;
               h21 = zero;
               h22 = zero;
               d1 = zero;
               d2 = zero;
               x1 = zero;
            }
            else
            {
               flag = one;
               h11 = p1 / p2;
               h22 = x1 / y1;
               u = one + h11 * h22;
               temp = d2 / u;
               d2 = d1 / u;
               d1 = temp;
               x1 = y1 * u;
            }
         }
         if(d1 != zero)
         {
            while((d1<=rgamsq)||(d1>=gamsq))
            {
               if (flag == zero)
               {
                  h11 = one;
                  h22 = one;
                  flag = -one;
               }
               else
               {
                  h21 = -one;
                  h12 = one;
                  flag = -one;
               }
               if(d1<=rgamsq)
               {
                  d1 = d1*gam*gam;
                  x1 = x1/gam;
                  h11 = h11/gam;
                  h12 = h12/gam;
               }
               else
               {
                  d1 = d1/(gam*gam);
                  x1 = x1*gam;
                  h11 = h11*gam;
                  h12 = h12*gam;
               }
            }
         }
         if(d2 != zero)
         {
            while((abs(d2)<=rgamsq)||(abs(d2)>=gamsq))
            {
               if(flag == zero)
               {
                  h11 = one;
                  h22 = one;
                  flag = -one;
               }
               else
               {
                  h21 = -one;
                  h12 = one;
                  flag = -one;
               }
               if(abs(d2)<=rgamsq)
               {
                  d2 = d2*gam*gam;
                  h21 = h21/gam;
                  h22 = h22/gam;
               }
               else
               {
                  d2 = d2/(gam*gam);
                  h21 = h21*gam;
                  h22 = h22*gam;
               }
            }
         }
      }
      if(flag < zero)
      {
         param[1] = h11;
         param[2] = h21;
         param[3] = h12;
         param[4] = h22;
      }
      else if(flag == zero)
      {
         param[2] = h21;
         param[3] = h12;
      }
      else
      {
         param[1] = h11;
         param[4] = h22;
      }
      param[0]=flag;
      return 0;
   }
}
#endif
