//
//  larfg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/12/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _larfg_h
#define _larfg_h

/// @file larfg.h Generates an elementary Householder reflection.

#include <cmath>
#include "lamch.h"
#include "nrm2.h"
#include "lapy2.h"
#include "lapy3.h"
#include "scal.h"
#include "ladiv.h"
#include "latl.h"

namespace LATL
{
   /// @brief Generates a real elementary Householder reflection.
   ///
   /// The real version of larfg generates a real elementary Householder reflection H of order n, such that
   /// 
   ///        H * ( alpha ) = ( beta ),   H' * H = I.
   ///            (   x   )   (   0  )
   /// 
   /// where alpha and beta are scalars, and x is an (n-1)-element vector. H is represented in the form
   /// 
   ///        H = I - tau * ( 1 ) * ( 1 v' ) 
   ///                      ( v )
   /// 
   /// where tau is a real scalar and v is a real (n-1)-element vector.
   /// Note that H is symmetric.
   /// 
   /// If the elements of x are all zero, then tau = 0 and H is taken to be
   /// the identity matrix.
   /// 
   /// Otherwise  1 <= tau <= 2.
   ///
   /// @tparam real_t Floating point type.
   /// @param n The order of the elementary Householder reflection.
   /// @param alpha On entry, the value alpha.  On exit, it is overwritten with the value beta.
   /// @param x Array of length 1+(n-2)*abs(incx).  On entry, the vector x.  On exit, it is overwritten with the vector v.
   /// @param incx  The increment between elements of x; incx > 0.
   /// @param tau On exit, the value tau.
   /// @ingroup AUX
   
   template<typename real_t>
   void LARFG(int_t n, real_t &alpha, real_t *x, int_t incx, real_t &tau)
   {
      using std::abs;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t safemin=LAMCH<real_t>('S')/LAMCH<real_t>('E');
      const real_t rsafemin=one/safemin;
      
      tau=zero;
      if(n>0)
      {
         int_t knt=0;
         real_t xnorm=NRM2<real_t>(n-1,x,incx);
         if(xnorm>zero)
         {
            real_t temp=abs(LAPY2(alpha,xnorm));
            real_t beta=(alpha<zero)?temp:-temp;
            if(abs(beta)<safemin)
            {
               while(abs(beta)<safemin)
               {
                  knt++;
                  SCAL<real_t>(n-1,rsafemin,x,incx);
                  beta*=rsafemin;
                  alpha*=rsafemin;
               }
               xnorm=NRM2<real_t>(n-1,x,incx);
               temp=abs(LAPY2(alpha,xnorm));
               beta=(alpha<zero)?temp:-temp;
            }
            tau=(beta-alpha)/beta;
            SCAL<real_t>(n-1,one/(alpha-beta),x,incx);
            for(int_t j=0;j<knt;j++)
               beta*=safemin;
            alpha=beta;
         }
      }
   }
   
   /// @brief Generates a complex elementary Householder reflector.
   ///
   /// The complex version of larfg generates a complex elementary Householder reflector H of order n, such that
   /// 
   ///        H' * ( alpha ) = ( beta ),   H' * H = I.
   ///             (   x   )   (   0  )
   /// 
   /// where alpha and beta are scalars, with beta real, and x is an (n-1)-element complex vector. 
   /// H is represented in the form
   /// 
   ///        H = I - tau * ( 1 ) * ( 1 v' ) 
   ///                      ( v )
   /// 
   /// where tau is a complex scalar and v is a complex (n-1)-element vector.
   /// Note that H is not hermitian.
   /// 
   /// If the elements of x are all zero and alpha is real, then tau = 0 and 
   /// H is taken to be the identity matrix.
   /// 
   /// Otherwise  1 <= real(tau) <= 2 and abs(tau-1) <= 1.
   ///
   /// @tparam real_t Floating point type.
   /// @param n The order of the elementary Householder reflection.
   /// @param alpha On entry, the value alpha.  On exit, it is overwritten with the value beta.
   /// @param x Array of length 1+(n-2)*abs(incx).  On entry, the vector x.  On exit, it is overwritten with the vector v.
   /// @param incx  The increment between elements of x; incx > 0.
   /// @param tau On exit, the value tau.
   /// @ingroup AUX
   
   template<typename real_t>
   void LARFG(int_t n, complex<real_t> &alpha, complex<real_t> *x, int_t incx, complex<real_t> &tau)
   {
      using std::real;
      using std::imag;
      using std::abs;

      const real_t one(1.0);
      const real_t zero(0.0);
      const complex<real_t> czero(0.0,0.0);
      const real_t safemin=LAMCH<real_t>('S')/LAMCH<real_t>('E');
      const real_t rsafemin=one/safemin;
      
      tau=czero;
      if(n>0)
      {
         int_t knt=0;      
         real_t xnorm=NRM2<real_t>(n-1,x,incx);
         if((xnorm!=zero)||(imag(alpha)!=zero))
         {
            real_t temp=abs(LAPY3(real(alpha),imag(alpha),xnorm));
            real_t beta=(real(alpha)<zero)?temp:-temp;
            if(abs(beta)<safemin)
            {
               while(abs(beta)<safemin)
               {
                  knt++;
                  SCAL<real_t>(n-1,rsafemin,x,incx);
                  beta*=rsafemin;
                  alpha*=rsafemin;
               }
               xnorm=NRM2<real_t>(n-1,x,incx);
               temp=abs(LAPY3(real(alpha),imag(alpha),xnorm));
               beta=(real(alpha)<zero)?temp:-temp;
            }
            tau=complex<real_t>((beta-real(alpha))/beta,-imag(alpha)/beta);
            alpha=LADIV<real_t>(complex<real_t>(one,zero),alpha-beta);
            SCAL<real_t>(n-1,alpha,x,incx);
            for(int_t j=0;j<knt;j++)
               beta*=safemin;
            alpha=beta;
         }
      }
   }
}

#endif
