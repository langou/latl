//
//  larft.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 2/11/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _larft_h
#define _larft_h

/// @file larft.h Forms the triangular factor T of a block reflector.

#include <cctype>
#include <algorithm>
#include "gemv.h"
#include "gemm.h"
#include "trmv.h"
#include "latl.h"

namespace LATL
{
   /// @brief Forms the triangular factor T of a real block reflector H of order n,
   /// which is defined as a product of k elementary reflectors.
   ///
   ///               If direct = 'F', H = H_1 H_2 . . . H_k and T is upper triangular.
   ///               If direct = 'B', H = H_k . . . H_2 H_1 and T is lower triangular.
   ///
   ///  If storeV = 'C', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th column of the array V, and
   ///
   ///               H  =  I - V * T * V'
   ///
   ///  If storeV = 'R', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th row of the array V, and
   ///
   ///               H  =  I - V' * T * V
   ///
   ///  The shape of the matrix V and the storage of the vectors which define
   ///  the H(i) is best illustrated by the following example with n = 5 and
   ///  k = 3. The elements equal to 1 are not stored.
   ///
   ///               direct='F' & storeV='C'          direct='F' & storeV='R'
   ///               -----------------------          -----------------------
   ///               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
   ///                   ( v1  1    )                     (     1 v2 v2 v2 )
   ///                   ( v1 v2  1 )                     (        1 v3 v3 )
   ///                   ( v1 v2 v3 )
   ///                   ( v1 v2 v3 )
   ///
   ///               direct='B' & storeV='C'          direct='B' & storeV='R'
   ///               -----------------------          -----------------------
   ///               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
   ///                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
   ///                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
   ///                   (     1 v3 )
   ///                   (        1 )
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param direct Specifies the direction in which the elementary reflectors are multiplied to form the block reflector.
   ///
   ///               'F' : forward
   ///               'B' : backward
   ///
   /// @param storeV Specifies how the vectors which define the elementary reflectors are stored.
   ///
   ///               'C' : columnwise
   ///               'R' : rowwise
   ///
   /// @param n The order of the block reflector H. n >= 0.
   /// @param k The order of the triangular factor T, or the number of elementary reflectors. k >= 1.
   /// @param[in] V Real matrix containing the vectors defining the elementary reflector H.
   /// If stored columnwise, V is n-by-k.  If stored rowwise, V is k-by-n.
   /// @param ldV Column length of the matrix V.  If stored columnwise, ldV >= n.
   /// If stored rowwise, ldV >= k.
   /// @param[in] tau Real vector of length k containing the scalar factors of the elementary reflectors H.
   /// @param[out] T Real matrix of size k-by-k containing the triangular factor of the block reflector.
   /// If the direction of the elementary reflectors is forward, T is upper triangular;
   /// if the direction of the elementary reflectors is backward, T is lower triangular.
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
                  GEMV<real_t>('T',n-i-1,i,-tau[i],&V0[i+1],ldV,&V[i+1],1,one,T,1);
               }
               else // storeV=='R'
               {
                  for(int_t j=0;j<i;j++)
                     T[j]=-tau[i]*V[j];
                  GEMV<real_t>('N',i,n-i-1,-tau[i],V+ldV,ldV,V+i+ldV,ldV,one,T,1);
               }
               TRMV<real_t>('U','N','N',i,T0,ldT,T,1);
               T[i]=tau[i];
            }
            T+=ldT;
            V+=ldV;
         }
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
                  GEMV<real_t>('T',n-k+i,k-i-1,-tau[i],V+ldV,ldV,V,1,one,T+i+1,1);
               }
               else // storeV=='R'
               {
                  real_t *v=V+(n-k)*ldV;
                  for(int_t j=i+1;j<k;j++)
                  {
                     T[j]=-tau[i]*v[j];
                  }
                  GEMV<real_t>('N',k-i-1,n-k+i,-tau[i],V0+i+1,ldV,V0+i,ldV,one,T+i+1,1);
               }
               TRMV<real_t>('L','N','N',k-i-1,T+i+1+ldT,ldT,T+i+1,1);
            }
            T[i]=tau[i];
         }
      }
      return 0;
   }

   /// @brief Forms the triangular factor T of a complex block reflector H of order n,
   /// which is defined as a product of k elementary reflectors.
   ///
   ///               If direct = 'F', H = H_1 H_2 . . . H_k and T is upper triangular.
   ///               If direct = 'B', H = H_k . . . H_2 H_1 and T is lower triangular.
   ///
   ///  If storeV = 'C', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th column of the array V, and
   ///
   ///               H  =  I - V * T * V'
   ///
   ///  If storeV = 'R', the vector which defines the elementary reflector
   ///  H(i) is stored in the i-th row of the array V, and
   ///
   ///               H  =  I - V' * T * V
   ///
   ///  The shape of the matrix V and the storage of the vectors which define
   ///  the H(i) is best illustrated by the following example with n = 5 and
   ///  k = 3. The elements equal to 1 are not stored.
   ///
   ///               direct='F' & storeV='C'          direct='F' & storeV='R'
   ///               -----------------------          -----------------------
   ///               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
   ///                   ( v1  1    )                     (     1 v2 v2 v2 )
   ///                   ( v1 v2  1 )                     (        1 v3 v3 )
   ///                   ( v1 v2 v3 )
   ///                   ( v1 v2 v3 )
   ///
   ///               direct='B' & storeV='C'          direct='B' & storeV='R'
   ///               -----------------------          -----------------------
   ///               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
   ///                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
   ///                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
   ///                   (     1 v3 )
   ///                   (        1 )
   ///
   /// @param direct Specifies the direction in which the elementary reflectors are multiplied to form the block reflector.
   ///
   ///               'F' : forward
   ///               'B' : backward
   ///
   /// @param storeV Specifies how the vectors which define the elementary reflectors are stored.
   ///
   ///               'C' : columnwise
   ///               'R' : rowwise
   ///
   /// @param n The order of the block reflector H. n >= 0.
   /// @param k The order of the triangular factor T, or the number of elementary reflectors. k >= 1.
   /// @param[in] V Complex matrix containing the vectors defining the elementary reflector H.
   /// If stored columnwise, V is n-by-k.  If stored rowwise, V is k-by-n.
   /// @param ldV Column length of the matrix V.  If stored columnwise, ldV >= n.
   /// If stored rowwise, ldV >= k.
   /// @param[in] tau Complex vector of length k containing the scalar factors of the elementary reflectors H.
   /// @param[out] T Complex matrix of size k-by-k containing the triangular factor of the block reflector.
   /// If the direction of the elementary reflectors is forward, T is upper triangular;
   /// if the direction of the elementary reflectors is backward, T is lower triangular.
   /// @param ldT Column length of the matrix T.  ldT >= k.
   
   template<typename real_t>
   int larft(char direct, char storeV, int_t n, int_t k, complex<real_t> *V, int_t ldV, complex<real_t> *tau, complex<real_t> *T, int_t ldT)
   {
      const complex<real_t> one(1.0);
      const complex<real_t> zero(0.0);
      using std::toupper;
      using std::max;
      using std::min;
      using std::conj;
      
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
      
      complex<real_t> *V0=V;
      complex<real_t> *T0=T;
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
                  complex<real_t> *v=V0;
                  for(int_t j=0;j<i;j++)
                  {
                     T[j]=-tau[i]*conj(v[i]);
                     v+=ldV;
                  }
                  GEMV<real_t>('C',n-i-1,i,-tau[i],&V0[i+1],ldV,&V[i+1],1,one,T,1);
               }
               else // storeV=='R'
               {
                  for(int_t j=0;j<i;j++)
                     T[j]=-tau[i]*V[j];
                  GEMM<real_t>('N','C',i,1,n-i-1,-tau[i],V+ldV,ldV,V+i+ldV,ldV,one,T,ldT);
               }
               TRMV<real_t>('U','N','N',i,T0,ldT,T,1);
               T[i]=tau[i];
            }
            T+=ldT;
            V+=ldV;
         }
      }
      else // direct=='B'
      {
         complex<real_t> *T0=T;
         complex<real_t> *V0=V;
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
                  complex<real_t> *v=V+ldV;
                  for(int_t j=i+1;j<k;j++)
                  {
                     T[j]=-tau[i]*conj(v[n-k+i]);
                     v+=ldV;
                  }
                  GEMV<real_t>('C',n-k+i,k-i-1,-tau[i],V+ldV,ldV,V,1,one,T+i+1,1);
               }
               else // storeV=='R'
               {
                  complex<real_t> *v=V+(n-k)*ldV;
                  for(int_t j=i+1;j<k;j++)
                  {
                     T[j]=-tau[i]*v[j];
                  }
                  GEMM<real_t>('N','C',k-i-1,1,n-k+i,-tau[i],V0+i+1,ldV,V0+i,ldV,one,T+i+1,ldT);
               }
               TRMV<real_t>('L','N','N',k-i-1,T+i+1+ldT,ldT,T+i+1,1);
            }
            T[i]=tau[i];
         }
      }
      return 0;
   }
}

#endif
