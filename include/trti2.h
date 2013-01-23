//
//  triti2.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _triti2_h
#define _triti2_h

#include "latl.h"
#include "trmv.h"
#include "scal.h"

namespace latl
{
   template <typename real_t>
   int_t trti2(char uplo, char diag, int_t n, real_t *A, int_t ldA)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const real_t one=1.0;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      if((diag!='N')&&(diag!='U'))
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      else if(n==0)
         return 0;
      
      real_t AJJ;
      if(uplo=='U')
      {
         for(int_t j=0;j<n;j++)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            int_t info=trmv<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            if(info!=0) return info;
            scal<real_t>(j,AJJ,A+j*ldA,1);
         }
      }
      else
      {
         for(int_t j=n-1;j>=0;j--)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            if(j<n-1)
            {
               int_t info=trmv('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               if(info!=0) return info;
               scal<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }

   template <typename real_t>
   int_t trti2(char uplo, char diag, int_t n, complex<real_t> *A, int_t ldA)
   {
      using std::toupper;
      using std::real;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const complex<real_t> one=1.0;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      if((diag!='N')&&(diag!='U'))
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      else if(n==0)
         return 0;
      
      complex<real_t> AJJ;
      if(uplo=='U')
      {
         for(int_t j=0;j<n;j++)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            int_t info=trmv<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            if(info!=0) return info;
            scal<real_t>(j,AJJ,A+j*ldA,1);
         }
      }
      else
      {
         for(int_t j=n-1;j>=0;j--)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            if(j<n-1)
            {
               int_t info=trmv('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               if(info!=0) return info;
               scal<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }
}

#endif
