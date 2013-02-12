//
//  larft.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 2/11/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _larft_h
#define _larft_h

/// @file larft.h Forms the triangular factor T of a block reflector H = I - vtvH.

#include <cctype>
#include <algorithm>
#include "gemv.h"
#include "trmv.h"
#include "latl.h"

namespace latl
{
   /// @brief Forms the triangular factor T of a real block reflector H of order n,
   /// which is defined as a product of k elementary reflectors.
   ///
   ///  If direct = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
   ///  If direct = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
   ///  If storeV = 'C', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th column of the array V, and
   ///
   ///         H  =  I - V * T * V'
   ///
   ///  If storeV = 'R', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th row of the array V, and
   ///
   ///         H  =  I - V' * T * V
   ///
   ///
   ///  The shape of the matrix V and the storage of the vectors which define
   ///  the H(i) is best illustrated by the following example with n = 5 and
   ///  k = 3. The elements equal to 1 are not stored.
   ///
   ///  direct = 'F' and storeV = 'C':         direct = 'F' and storeV = 'R':
   ///
   ///               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
   ///                   ( v1  1    )                     (     1 v2 v2 v2 )
   ///                   ( v1 v2  1 )                     (        1 v3 v3 )
   ///                   ( v1 v2 v3 )
   ///                   ( v1 v2 v3 )
   ///
   ///  direct = 'B' and storeV = 'C':         direct = 'B' and storeV = 'R':
   ///
   ///               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
   ///                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
   ///                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
   ///                   (     1 v3 )
   ///                   (        1 )
   ///
   /// @param direct Specifies the order in which the elementary reflectors are
   ///               multiplied to form the block reflector:
   ///
   ///               direct = 'F': H = H(1) H(2) . . . H(k) (Forward)
   ///               direct = 'B': H = H(k) . . . H(2) H(1) (Backward)
   ///
   /// @param storeV Specifies how the vectors which define the elementary
   ///               reflectors are stored.
   ///
   ///               storeV = 'C': columnwise
   ///               storeV = 'R': rowwise
   ///
   /// @param n The order of the block reflector H. n >= 0.
   /// @param k The order of the triangular factor T (= the number of elementary reflectors). k >= 1.
   /// @param V Real matrix containing the vectors defining the elementary reflector H.
   ///
   ///               storeV = 'C' : V is n-by-k
   ///               storeV = 'R' : V is k-by-n
   ///
   /// @param ldV Column length of the matrix V.
   ///
   ///               storeV = 'C' : ldV >= n
   ///               storeV = 'R' : ldV >= k
   ///
   /// @param tau Real vector of length k containing the scalar factors of the elementary reflectors H.
   ///
   /// @param T Real matrix of size k-by-k containing the triangular factor of the block reflector.
   ///
   ///               direct = 'F' :  T is upper triangular
   ///               direct = 'B' :  T is lower triangular.
   ///
   /// @param ldT Column length of the matrix T.  ldT >= k.

   template<typename real_t>
   int larft(char direct, char storeV, int_t n, int_t k, real_t *V, int_t ldV, real_t *tau, real_t *T, int_t ldT)
   {
      const real_t one(1.0);
      const real_t zero(0.0);
      using std::toupper;
      using std::max;
      using std::min;

      direct=toupper(direct);
      storeV=toupper(storeV);

      if((direct!='F')&&(direct!='B'))
         return -1;
      else if((storeV!='C')&&(storeV!='R'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<1)
         return -3;
      else if(((storeV=='C')&&(ldV<n))||((storeV=='R')&&(ldV<k)))
         return -5;
      else if(ldT<k)
         return -9;

      if(n==0)
         return 0;

      real_t *V0=V;
      real_t *T0=T;
      if(direct=='F')
      {
         for(int_t i=0;i<k;i++)
         {
            if(tau[i]==zero)
            {
               for(int j=0;j<=i;j++)
                  T[j]=zero;
            }
            else
            {
               if(storeV=='C')
               {
                  real_t *v=V0;
                  for(int_t j=0;j<i;j++)
                  {
                     T[j]=-tau[i]*v[i];
                     v+=ldV;
                  }
                  gemv<real_t>('T',n-i,i,-tau[i],V0+i+1,ldV,V+i+1,1,one,T,1);

               }
               else // storeV=='R'
               {
                  for(int_t j=0;j<i;j++)
                     T[j]=-tau[i]*V[j];
                  gemv<real_t>('N',i,n-i,-tau[i],V+ldV,ldV,V+i+ldV,ldV,one,T,1);
               }
               trmv<real_t>('U','N','N',i,T0,ldT,T,1);
               T[i]=tau[i];
            }
         }
         T+=ldT;
         V+=ldV;
      }
      else // direct=='B'
      {
         real_t *T0=T;
         real_t *V0=V;
         T+=(k-1)*ldT;
         V+=(k-1)*ldV;
         T[k-1]=tau[k-1];
         for(int_t i=k-2;i>=0;--i)
         {
            T-=ldT;
            V-=ldV;
            if(tau[i]==zero)
            {
               for(int_t j=i;j<k;j++)
                  T[j]=zero;
            }
            else
            {
               if(storeV=='C')
               {
                  real_t *v=V+ldV;
                  for(int_t j=i+1;j<k;j++)
                  {
                     T[j]=-tau[i]*v[n-k+i];
                     v+=ldV;
                  }
                  gemv<real_t>('T',n-k+i,k-i,-tau[i],V+ldV,ldV,V,1,one,T+i+1,1);
               }
               else // storeV=='R'
               {
                  real_t *v=V+(n-k+1)*ldV;
                  for(int_t j=i+1;j<k;j++)
                  {
                     T[j]=-tau[i]*v[j];
                     v+=ldV;
                  }
                  gemv<real_t>('N',k-i,n-k+i,-tau[i],V0+i+1,ldV,V+i,ldV,one,T+i+1,1);
               }
               trmv<real_t>('L','N','N',k-i,T+i+1+ldT,ldT,T+i+1,1);
            }
            T[i]=tau[i];
         }
      }
      return 0;
   }
}

#endif
