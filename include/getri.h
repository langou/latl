//
//  getri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/23/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _getri_h
#define _getri_h

#include "latl.h"
#include "trtri.h"
#include "gemv.h"
#include "gemm.h"
#include "trsm.h"

namespace latl
{
   template <typename real_t>
      int_t getri(int_t n, real_t *A, int_t ldA, int_t *ipiv, int_t nb=32)
      {
         const real_t zero=0.0;
         const real_t one=1.0;

         if(n<0)
            return -1;
         else if(ldA<n)
            return -3;
         else if(n==0)
            return 0;

         real_t *work=new real_t[n*nb];

         int_t info=trtri<real_t>('U','N',n,A,ldA);
         if(info>0)
            return info;

         if(nb>=n)
         {
            for(int_t j=n-1;j>=0;j--)
            {
               for(int_t i=j+1;i<n;i++)
               {
                  work[i]=A[i+j*ldA];
                  A[i+j*ldA]=zero;
               }
               if(j<n-1)
                  gemv<real_t>('N', n, n-j-1, -one, A+(j+1)*ldA, ldA, work+j+1, 1, one, A+j*ldA, 1);
            }
         }
         else
         {
            int_t nn = ( (n-1) / nb ) * nb + 1;  // LINE 223 OF DGETRI.F
            for(int_t j=nn;j>=1;j-=nb)
            {
               int_t jb = std::min( nb, n-j+1 );
               for(int_t jj=j; jj<=j+jb-1; jj++)
               {
                  for(int_t i=jj+1;i<=n;i++)
                  {
                     work[(i-1)+(jj-j)*n] = A[(i-1)+(jj-1)*ldA];
                     A[(i-1)+(jj-1)*ldA] = zero;
                  }
               }
 
               if( j+jb <= n)
                  gemm<real_t>('N', 'N', n, jb, n-j-jb+1, -one, A+(j+jb-1)*ldA, ldA, work+j+jb-1, n, one, A+(j-1)*ldA, ldA);

               trsm<real_t>('R', 'L', 'N', 'U', n, jb, one, work+(j-1), n, A+(j-1)*ldA, ldA);
            }
         }

         for(int_t j=n-2;j>=0;j--)
         {
            int_t jp = ipiv[j];
            if(jp != j)
               swap(n,A+j*ldA,1,A+jp*ldA,1);
         }

         delete [] work;
         return 0;
      }
}

#endif

