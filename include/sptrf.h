//
//  sptrf.h
//  LAPACK Template Library
//
//  Created by Stephanie Patterson on 8/2/12.
//
//

#ifndef _sptrf_h
#define _sptrf_h

/// @file sptrf.h Computes the factorization of a symmetric matrix stored in packed format.

#include "spr.h"
#include "swap.h"
#include "scal.h"
#include "imax.h"
#include "latl.h"
#include <iostream>

namespace LATL
{
   template<typename real_t>
   int_t sptrf( const int_t n, real_t * const AP, int_t * const IPIV, bool * const BSDV)
   {
      /// @brief Computes the factorization of a real symmetric matrix A stored in packed format using the Bunch-Kaufman diagonal pivoting method:
      ///
      ///     A = U * D * U'  
      ///
      /// where U is a product of permutation and unit upper triangular matrices, and D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
      /// @return 0 if success.
      /// @return -i if the ith argument is invalid.
      /// @return i+1 if the ith diagonal element is exactly zero.  The factorization has been completed, but the block diagonal matrix D is singular and division by zero will occur if it is used to solve a system of equations.
      /// @tparam real_t Real floating point type.
      /// @param n Order of the matrix A.  n >= 0
      /// @param AP Real array size n-by-(n+1)/2.  On entry, the upper triangle of the symmetric matrix A, packed columnwise in a linear array.  The jth column begins at AP + sum(0, j) = AP + j * (j+1)/2.
      /// On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L, stored as a packed triangular matrix.
      /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges of D.
      /// @param BSDV Bool array size n. On exit, contains the details of the block structure of D.  If BSDV[k] = 0, then rows and columns k and IPIV[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If BSDV[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if IPIV[k] = IPIV[k-1], then rows and columns k-1 and IPIV[k] were interchanged.
      /// @ingroup TRF
      
      if (n < 0)
         return -1;
      
      if ( n == 0)
         return 0;
      
      using std::isnan;
      using std::abs;
      using std::max;
      using std::sqrt;
      const real_t alpha = (1.0 + sqrt(17.0))/8.0;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t absakk, colmax, d11, d22, d12, r1, rowmax, temp, *APimax, *APk, *APkp, *APt, *APkm1, wkm1, wk;
      int_t imax, jmax, kstep, kp, kk, info = 0;
      
      int_t k = n-1;
         
      while (k >= 0)
      {
         kstep = 1;
         APk = AP + k*(k+1)/2;
         APkm1 = APk - k;
         
         absakk = abs(APk[k]);
         
         if (k > 0)
         {
            imax = LATL::imax(k, APk, 1);
            colmax = abs(APk[imax]);
         }
         else
         {
            colmax = zero;
         }
         
         if (max(absakk, colmax) == zero || isnan(absakk))
         {
            if (info == 0)
               info = k+1;
            kp = k;
         }
         else
         {
            if (absakk >= alpha*colmax)
            {
               kp = k;
            }
            else
            {
               rowmax = zero;
               jmax = imax;
               APimax = AP + imax*(imax+1)/2;
               for (int_t j = imax+1; j <= k; ++j)
               {
                  APimax += j;
                  if (abs(APimax[imax]) > rowmax)
                  {
                     rowmax = abs(APimax[imax]);
                     jmax = j;
                  }
               }
               APimax = AP + imax*(imax+1)/2;
               if (imax > 0)
               {
                  jmax = LATL::imax(imax, APimax, 1);
                  rowmax = max(rowmax, abs(APimax[jmax]));
               }
                  
               if (absakk >= alpha*colmax*(colmax/rowmax))
               {
                  kp = k;
               }
               else if (abs(APimax[imax]) >= alpha*rowmax)
               {
                  kp = imax;
               }
               else
               {
                  kp = imax;
                  kstep = 2;
               }
            }
               
            kk = k - kstep + 1;
            if (kstep == 2)
            {
               APkp = AP + k*(k-1)/2;
            }
            else
            {
               APkp = APk;
            }
            
            if (kp != kk)
            {
               LATL::swap(kp, APkp, 1, APimax, 1);
               APt = APimax;
               for (int_t j = kp+1; j < kk; ++j)
               {
                  APt += j;
                  temp = APkp[j];
                  APkp[j] = APt[kp];
                  APt[kp] = temp;
               }
               temp = APkp[kk];
               APkp[kk] = APimax[kp];
               APimax[kp] = temp;
               if (kstep == 2)
               {
                  temp = APk[k-1];
                  APk[k-1] = APk[kp];
                  APk[kp] = temp;
               }
            }
            
            if (kstep == 1)
            {
               r1 = one/ APk[k];
               LATL::spr('U', k, -r1, APk, 1, AP);
               LATL::scal(k, r1, APk, 1);
            }
            else
            {
               if (k > 1)
               {
                  d12 = APk[k-1];
                  d22 = APkm1[k-1]/d12;
                  d11 = APk[k]/d12;
                  d12 = (one/(d11*d22-one))/d12;
                  
                  for (int_t j = k-2; j >= 0; --j)
                  {
                     wkm1 = d12*(d11*APkm1[j] - APk[j]);
                     wk = d12*(d22*APk[j] - APkm1[j]);
                     APt = AP + j*(j+1)/2;
                     for (int_t i = j; i >= 0; --i)
                     {
                        APt[i] -= (APk[i]*wk + APkm1[i]*wkm1);
                     }
                     APk[j] = wk;
                     APkm1[j] = wkm1;
                  }
               }
            }
         }
         if (kstep == 1)
         {
            IPIV[k] = kp;
            BSDV[k] = 0;
         }
         else
         {
            IPIV[k] = kp;
            IPIV[k-1] = kp;
            BSDV[k] = 1;
            BSDV[k-1] = 1;
         }
         
         k -= kstep;
      }
      return info;
   }
   
   template<typename real_t>
   int_t sptrf( const int_t n, complex<real_t> * const AP, int_t * const IPIV, bool * const BSDV)
   {
      /// @brief Computes the factorization of a complex symmetric matrix A stored in packed format using the Bunch-Kaufman diagonal pivoting method:
      ///
      ///     A = U * D * U'
      ///
      /// where U is a product of permutation and unit upper triangular matrices, and D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
      /// @return 0 if success.
      /// @return -i if the ith argument is invalid.
      /// @return i+1 if the ith diagonal element is exactly zero.  The factorization has been completed, but the block diagonal matrix D is singular and division by zero will occur if it is used to solve a system of equations.
      /// @tparam real_t Real floating point type.
      /// @param n Order of the matrix A.  n >= 0
      /// @param AP Complex array size n-by-(n+1)/2.  On entry, the upper triangle of the symmetric matrix A, packed columnwise in a linear array.  The jth column begins at AP + sum(0, j) = AP + j * (j+1)/2.
      /// On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L, stored as a packed triangular matrix.
      /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges of D.
      /// @param BSDV Bool array size n. On exit, contains the details of the block structure of D.  If BSDV[k] = 0, then rows and columns k and IPIV[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If BSDV[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if IPIV[k] = IPIV[k-1], then rows and columns k-1 and IPIV[k] were interchanged.
      /// @ingroup TRF
      
      if (n < 0)
         return -1;
      
      if ( n == 0)
         return 0;
      
      using std::isnan;
      using std::abs;
      using std::max;
      using std::sqrt;
      const real_t alpha = (1.0 + sqrt(17.0))/8.0;
      const real_t zero(0.0);
      const complex<real_t> one(1.0);
      real_t absakk, colmax, rowmax;
      complex<real_t> d11, d22, d12, r1, temp, *APimax, *APk, *APkp, *APt, *APkm1, wkm1, wk;
      int_t imax, jmax, kstep, kp, kk, info = 0;
      
      int_t k = n-1;
      
      while (k >= 0)
      {
         kstep = 1;
         APk = AP + k*(k+1)/2;
         APkm1 = APk - k;
         
         absakk = abs(real(APk[k]))+abs(imag(APk[k]));
         
         if (k > 0)
         {
            imax = LATL::imax(k, APk, 1);
            colmax = abs(real(APk[imax]))+abs(imag(APk[imax]));
         }
         else
         {
            colmax = zero;
         }
         
         if (max(absakk, colmax) == zero || isnan(absakk))
         {
            if (info == 0)
               info = k+1;
            kp = k;
         }
         else
         {
            if (absakk >= alpha*colmax)
            {
               kp = k;
            }
            else
            {
               rowmax = zero;
               jmax = imax;
               APimax = AP + imax*(imax+1)/2;
               for (int_t j = imax+1; j <= k; ++j)
               {
                  APimax += j;
                  if (abs(real(APimax[imax]))+abs(imag(APimax[imax])) > rowmax)
                  {
                     rowmax = abs(real(APimax[imax]))+abs(imag(APimax[imax]));
                     jmax = j;
                  }
               }
               APimax = AP + imax*(imax+1)/2;
               if (imax > 0)
               {
                  jmax = LATL::imax(imax, APimax, 1);
                  rowmax = max(rowmax, abs(real(APimax[jmax]))+abs(imag(APimax[jmax])));
               }
               
               if (absakk >= alpha*colmax*(colmax/rowmax))
               {
                  kp = k;
               }
               else if (abs(real(APimax[imax]))+abs(imag(APimax[imax])) >= alpha*rowmax)
               {
                  kp = imax;
               }
               else
               {
                  kp = imax;
                  kstep = 2;
               }
            }
            
            kk = k - kstep + 1;
            if (kstep == 2)
            {
               APkp = AP + k*(k-1)/2;
            }
            else
            {
               APkp = APk;
            }
            
            if (kp != kk)
            {
               LATL::swap(kp, APkp, 1, APimax, 1);
               APt = APimax;
               for (int_t j = kp+1; j < kk; ++j)
               {
                  APt += j;
                  temp = APkp[j];
                  APkp[j] = APt[kp];
                  APt[kp] = temp;
               }
               temp = APkp[kk];
               APkp[kk] = APimax[kp];
               APimax[kp] = temp;
               if (kstep == 2)
               {
                  temp = APk[k-1];
                  APk[k-1] = APk[kp];
                  APk[kp] = temp;
               }
            }
            
            if (kstep == 1)
            {
               r1 = one/ APk[k];
               LATL::spr('U', k, -r1, APk, 1, AP);
               LATL::scal(k, r1, APk, 1);
            }
            else
            {
               if (k > 1)
               {
                  d12 = APk[k-1];
                  d22 = APkm1[k-1]/d12;
                  d11 = APk[k]/d12;
                  d12 = (one/(d11*d22-one))/d12;
                  
                  for (int_t j = k-2; j >= 0; --j)
                  {
                     wkm1 = d12*(d11*APkm1[j] - APk[j]);
                     wk = d12*(d22*APk[j] - APkm1[j]);
                     APt = AP + j*(j+1)/2;
                     for (int_t i = j; i >= 0; --i)
                     {
                        APt[i] -= (APk[i]*wk + APkm1[i]*wkm1);
                     }
                     APk[j] = wk;
                     APkm1[j] = wkm1;
                  }
               }
            }
         }
         if (kstep == 1)
         {
            IPIV[k] = kp;
            BSDV[k] = 0;
         }
         else
         {
            IPIV[k] = kp;
            IPIV[k-1] = kp;
            BSDV[k] = 1;
            BSDV[k-1] = 1;
         }
         
         k -= kstep;
      }
      return info;
   }
}


#endif
