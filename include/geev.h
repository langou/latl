//
//  geev.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _geev_h
#define _geev_h

/// @file geev.h Solves a nonsymmetric eigenvalue problem.

#include "lartg.h"
#include "rot.h"
#include "imax.h"
#include "lapy2.h"
#include "nrm2.h"
#include "gehrd.h"
#include "hseqr.h"
#include "orghr.h"
#include "trecv.h"
#include "lacpy.h"
#include "latl.h"

namespace LATL
{
   template<typename real_t>
   int GEEV(bool leftV,bool rightV,int_t n,real_t *A,int_t ldA,real_t *wr,real_t *wi,real_t *VL,int_t ldVL,real_t *VR,int_t ldVR)
   {
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

      // reduce to upper Hessenberg form

      real_t *tau=new real_t[n];
      GEHRD(n,A,ldA,tau);

      // compute eigenvalues and eigenvectors

      if(leftV&&rightV)
      {
         LACPY('L',n,n,A,ldA,VL,ldVL);
         ORGHR(n,VL,ldVL,tau);
         info=HSEQR('S','V',n,A,ldA,wr,wi,VL,ldVL);
         LACPY('F',n,n,VL,ldVL,VR,ldVR);
         TREVC('B','B',n,A,ldA,VL,ldVL,VR,ldVR,n);
      }
      else if(leftV)
      {
         LACPY('L',n,n,A,ldA,VL,ldVL);
         ORGHR(n,VL,ldVL,tau);
         info=HSEQR('S','V',n,A,ldA,wr,wi,VL,ldVL);
         TREVC('L','B',n,A,ldA,VL,ldVL,VR,ldVR,n);
      }
      else if(rightV)
      {
         LACPY('L',n,n,A,ldA,VR,ldVR);
         ORGHR(n,VR,ldVR,tau);
         info=HSEQR('S','V',n,A,ldA,wr,wi,VR,ldVR);
         TREVC('R','B',n,A,ldA,VL,ldVL,VR,ldVR,n);
      }
      else // eigenvalues only
      {
         info=HSEQR('E','N',n,A,ldA,wr,wi,VR,ldVR);
      }

      // normalize eigenvectors while making largest real
      
      if(leftV)
      {
         real_t cs,sn,r;
         real_t *t=new real_t[n];
         real_t *V=VL;
         real_t *W=V+ldVL;
         for(int_t i=0;i<n;i++)
         {
            if(wi[i]==zero)
            {
               real_t scale=one/NRM2(n,V,1);
               SCAL(n,scale,V,1);
            }
            else if(wi[i]>zero)
            {
               real_t scale=one/LAPY2(NRM2(n,V,1),NRM2(n,W,1));
               SCAL(n,scale,V,1);
               SCAL(n,scale,W,1);
               for(int_t k=0;k<n;k++)
                  t[i]=V[k]*V[k]+W[k]*W[k];
               int_t k=IMAX(n,t,1);
               LARTG(V[k],W[k],cs,sn,r);
               ROT(n,V,1,W,1,cx,sn);
               W[k]=zero;
            }
            V+=ldVL;
            W+=ldVL;
         }
         delete [] t;
      }

      if(rightV)
      {
         real_t cs,sn,r;
         real_t *t=new real_t[n];
         real_t *V=VR;
         real_t *W=V+ldVR;
         for(int_t i=0;i<n;i++)
         {
            if(wi[i]==zero)
            {
               real_t scale=one/NRM2(n,V,1);
               SCAL(n,scale,V,1);
            }
            else if(wi[i]>zero)
            {
               real_t scale=one/LAPY2(NRM2(n,V,1),NRM2(n,W,1));
               SCAL(n,scale,V,1);
               SCAL(n,scale,W,1);
               for(int_t k=0;k<n;k++)
                  t[i]=V[k]*V[k]+W[k]*W[k];
               int_t k=IMAX(n,t,1);
               LARTG(V[k],W[k],cs,sn,r);
               ROT(n,V,1,W,1,cx,sn);
               W[k]=zero;
            }
            V+=ldVR;
            W+=ldVR;
         }
         delete [] t;
      }

      delete [] tau;
   }
}
#endif
