//
//  rotm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _rotm_h
#define _rotm_h

/// @file rotm.h Applies modified Givens rotation matix to a pair of vectors.

#include "latl.h"

namespace LATL
{
   /// @brief Applies a modified Givens rottaion matrix H to a pair of real vectors (x,y).
   ///
   /// The Givens rotation matrix H with flag obtained by LATL::rotmg are applied to a pair of vectors (x,y).
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x,y.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of vector y.
   /// @param flag Indicates the form of the Givens rotation matrix H, obtained from LATL::rotmg.
   /// @param H Pointer to 2-by-2 Givens rotation matrix obtained from LATL::rotmg.
   /// @ingroup ROT

   template <typename real_t>
   void rotm(int_t n, real_t *x, int_t incx, real_t *y, int_t incy, real_t flag, real_t *H)
   {
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t h11(H[0]);
      const real_t h21(H[1]);
      const real_t h12(H[2]);
      const real_t h22(H[3]);
      real_t w, z;
      int_t i,ix,iy;
      
      if(n>0) 
      {
         if((incx==incy)&&(incx==1))
         {
            if(flag==-one)
            {
               for(i=0;i<n;i++)
               {
                  w=x[i];
                  z=y[i];
                  x[i]=w*h11+z*h12;
                  y[i]=w*h21+z*h22;
               }
            }
            else if(flag==zero)
            {
               for(i=0;i<n;i++)
               {
                  w=x[i];
                  z=y[i];
                  x[i]=w+z*h12;
                  y[i]=w*h21+z;
               }
            }
            else if(flag==one)     
            {
               for(i=0;i<n;i++)
               {
                  w=x[i];
                  z=y[i];
                  x[i]=w*h11+z;
                  y[i]=-w+z*h22;
               }
            }
         }
         else
         {
            ix=(incx>0)?0:(1-n)*incx;
            iy=(incy>0)?0:(1-n)*incy;
            if(flag==-one)
            {
               for(i=0;i<n;i++)
               {
                  w=x[ix];
                  z=y[iy];
                  x[ix]=w*h11+z*h12;
                  y[iy]=w*h21+z*h22;
                  ix+=incx;
                  iy+=incy;
               }
            }
            else if(flag==zero)     
            {
               for(i=0;i<n;i++)
               {
                  w=x[ix];
                  z=y[iy];
                  x[ix]=w+z*h12;
                  y[iy]=w*h21+z;
                  ix+=incx;
                  iy+=incy;
               }
            }
            else if(flag==one)    
            {
               for(i=0;i<n;i++)
               {
                  w=x[ix];
                  z=y[iy];
                  x[ix]=w*h11+z;
                  y[iy]=-w+z*h22;
                  ix+=incx;
                  iy+=incy;
               }
            }
         }
      }
   }
}
#endif
   
