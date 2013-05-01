//
//  geev.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _geev_h
#define _geev_h

/// @file geev.h Solves nonsymmetric eigenvalue problem.

#include <cmath>
#include "gehrd.h"
#include "lange.h"
#include "lamch.h"
#include "labad.h"
#include "lascl.h"
#include "latl.h"

namespace LATL
{
   template<typename real_t>
   int GEEV(bool leftV,bool rightV,int_t n,real_t *A,int_t ldA,real_t *wr,real_t *wi,real_t *VL,int_t ldVL,real_t *VR,int_t ldVR)
   {
      using std::sqrt;
      const real_t one(1.0);
      const real_t zero(1.0);
      if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      else if(ldVL<n)
         return -9;
      else if(ldVR<n)
         return -11;
      else if(n==0)
         return 0;

      // machine constants

      const real_t eps=LAMCH('P');
      real_t small=LAMCH('S');
      real_t big=one/small;
      LABAD(small,big);
      small=sqrt(small)/eps;
      big=one/small;

      // scale A if necessary
      
      real_t Anorm=LANGE('M',n,n,A,ldA);
      real_t scale=one;
      if((Anorm>zero)&&(Anorm<small))
      {
         scale=small;
         LASCL('G',0,0,Anorm,scale,n,n,A,ldA);
      }
      else if(Anorm>big)
      {
         scale=big;
         LASCL('G',0,0,Anorm,scale,n,n,A,ldA);
      }

      // reduce to upper Hessenberg form

      real_t *tau=new real_t[n];
      GEHRD(n,A,ldA,tau);

      // compute eigenvectors

      


   }
}
#endif
