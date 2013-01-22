//
//  sytf2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/29/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sytf2_h
#define _sytf2_h

/// @file sytf2.h Computes the factorization of a symmetric matrix A.

#include "imax.h"
#include "syr.h"
#include "scal.h"
#include "swap.h"
#include <cmath>
#include "latl.h"

namespace latl
{
   /// @brief Computes the factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method:
   ///
   ///         A = U*D*U' if uplo = 'U'
   ///         A = L'*D*L if uplo = 'L'
   ///
   /// where U (or L) is a product of permutation and unit upper (lower) triangular matrices, U' is the transpose of U, and D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is exactly zero.  The factorization has been completed, but the block diagonal matrix D is exactly singular and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the symmetric matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA Column length of the array A.  ldA >= n
   /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges of D.
   /// @param BSDV Bool array size n.  On exit, contains the details of the block structure of D.  If BSDV[k] = 0, then rows and columns k and IPIV[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If BSDV[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if uplo = 'U', and IPIV[k] = IPIV[k-1], then rows and columns k-1 and IPIV[k] were interchanged.  If uplo = 'L' and IPIV[k] = IPIV[k+1], then rows and columns k+1 and IPIV[k] were interchanged.
   
   template<typename real_t>
   int_t sytf2(const char uplo, const int_t n, real_t * const A, const int_t ldA, int_t * IPIV, bool * BSDV)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      const real_t alpha = (1.0 + std::sqrt(17.0))/8.0;
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t * Ak, absakk, colmax, * Akm1;
      int_t info = 0, imax = 0, kp;
      
      for (int_t i = 0; i < n; ++i)
      {
         BSDV[i] = 0;
      }
      
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kstep;
         while (k >= 0)
         {
            kstep = 1;
            Ak = A + ldA*k;
            Akm1 = Ak-ldA;
            
            absakk = std::abs(Ak[k]);
            
            if (k > 0)
            {
               imax = latl::imax(k, Ak, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == zero || std::isnan(absakk))
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
                  int_t jmax;
                  real_t rowmax;
                  
                  jmax = imax+1+ latl::imax(k-imax, A+ldA*(imax+1)+imax, ldA);
                  real_t * Ajmax = A + ldA*jmax;
                  real_t * Aimax = A + ldA*imax;
                  rowmax = std::abs(Ajmax[imax]);
                  if (imax > 0)
                  {
                     jmax = latl::imax(imax, Aimax, 1);
                     rowmax = std::max(rowmax, std::abs(Aimax[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (std::abs(Aimax[imax]) >= alpha*rowmax)
                     kp = imax;
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k-kstep+1;
               real_t * Akinc = A + ldA*kinc;
               real_t * Akp = A + ldA*kp;
               real_t r1;
               
               if (kp != kinc)
               {
                  latl::swap(kp, Akinc, 1, Akp, 1);
                  latl::swap(kinc-kp-1, Akinc+kp+1, 1, A+ldA*(kp+1)+kp, ldA);
                  real_t temp = Akinc[kinc];
                  Akinc[kinc] = Akp[kp];
                  Akp[kp] = temp;
                  if (kstep == 2)
                  {
                     temp = Ak[k-1];
                     Ak[k-1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               
               if (kstep == 1)
               {
                  r1 = one/Ak[k];
                  latl::syr(uplo, k, -r1, Ak, 1, A, ldA);
                  latl::scal(k, r1, Ak, 1);
               }
               else
               {
                  if (k > 1)
                  {
                     real_t d12, d22, d11;
                     d12 = Ak[k-1];
                     d22 = Akm1[k-1]/d12;
                     d11 = Ak[k]/d12;
                     d12 = (one/(d11*d22-one))/d12;
                     
                     real_t wk, wkm1;
                     for (int_t j = k-2; j >= 0; --j)
                     {
                        real_t * Aj = A+ldA*j;
                        wkm1 = d12*(d11*Akm1[j]-Ak[j]);
                        wk = d12*(d22*Ak[j]-Akm1[j]);
                        for (int_t i = j; i >= 0; --i)
                        {
                           Aj[i] -= (Ak[i]*wk+Akm1[i]*wkm1);
                        }
                        Ak[j] = wk;
                        Akm1[j] = wkm1;
                     }
                  }
               }
            }
            if (kstep == 1)
            {
               IPIV[k] = kp;
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
      }
      else
      {
         int_t k = 0, kstep = 1;
         while (k < n)
         {
            kstep = 1;
            Ak = A + ldA*k;
            
            absakk = std::abs(Ak[k]);
            
            if (k < n-1)
            {
               imax = k +1 + latl::imax(n-k-1, Ak+k+1, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == zero || std::isnan(absakk))
            {
               if (info == 0)
                  info = k;
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
                  kp = k;
               else
               {
                  int_t jmax = k + latl::imax(imax-k, Ak+imax, ldA);
                  real_t * Ajmax = A + ldA*jmax;
                  real_t *Aimax = A +ldA*imax;
                  real_t rowmax = std::abs(Ajmax[imax]);
                  if (imax < n-1)
                  {
                     jmax = imax +latl::imax(n-imax-1, Aimax+imax+1, 1);
                     rowmax = std::max(rowmax, std::abs(Aimax[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (std::abs(Aimax[imax]) >= alpha*rowmax)
                  {
                     kp = imax;
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k + kstep -1;
               real_t * Akinc = A+ldA*kinc;
               real_t * Akp = A+ldA*kp;
               if (kp != kinc)
               {
                  if (kp < n-1)
                     latl::swap(n-kp-1, Akinc+kp+1, 1, Akp+kp+1, 1);
                  latl::swap(kp-kinc-1, Akinc+kinc+1, 1, Akinc+ldA+kp, ldA);
                  real_t temp = Akinc[kinc];
                  Akinc[kinc] = Akp[kp];
                  Akp[kp] = temp;
                  if (kstep == 2)
                  {
                     temp = Ak[k+1];
                     Ak[k+1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               if (kstep == 1)
               {
                  if (k < n-1)
                  {
                     real_t r1 = one/Ak[k];
                     latl::syr(uplo, n-k-1, -r1, Ak+k+1, 1, Ak+ldA+k+1, ldA);
                     latl::scal(n-k-1, r1, Ak+k+1, 1);
                  }
               }
               else
               {
                  if (k < n-2)
                  {
                     real_t d21, d11, d22, temp, wk, wkp1;
                     real_t * Akp1 = A+ldA*(k+1);
                     d21 = Ak[k+1];
                     d11 = Akp1[k+1]/d21;
                     d22 = Ak[k]/d21;
                     temp = one/(d11*d22-one);
                     d21 = temp/d21;
                     
                     for (int_t j = k+2; j < n; ++j)
                     {
                        wk = d21*(d11*Ak[j]-Akp1[j]);
                        wkp1 = d21*(d22*Akp1[j]-Ak[j]);
                        real_t * Aj = A + ldA*j;
                        
                        for (int_t i = j; i< n; ++i)
                        {
                           Aj[i] -= Ak[i]*wk+Akp1[i]*wkp1;
                        }
                        Ak[j] = wk;
                        Akp1[j] = wkp1;
                     }
                  }
               }
            }
            if (kstep == 1)
               IPIV[k] = kp;
            else
            {
               IPIV[k] = kp;
               IPIV[k+1] = kp;
               BSDV[k] = 1;
               BSDV[k+1] = 1;
            }
            
            k += kstep;
         }
      }
      return info;
   }
   
   /// @brief Computes the factorization of a complex symmetric matrix A using the Bunch-Kaufman diagonal pivoting method:
   ///
   ///         A = U*D*U' if uplo = 'U'
   ///         A = L'*D*L if uplo = 'L'
   ///
   /// where U (or L) is a product of permutation and unit upper (lower) triangular matrices, U' is the transpose of U, and D is symmetric and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is exactly zero.  The factorization has been completed, but the block diagonal matrix D is exactly singular and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Complex array size ldA-by-n.  On entry, the symmetric matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA Column length of the array A.  ldA >= n
   /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges of D.
   /// @param BSDV Bool array size n.  On exit, contains the details of the block structure of D.  If BSDV[k] = 0, then rows and columns k and IPIV[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If BSDV[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if uplo = 'U', and IPIV[k] = IPIV[k-1], then rows and columns k-1 and IPIV[k] were interchanged.  If uplo = 'L' and IPIV[k] = IPIV[k+1], then rows and columns k+1 and IPIV[k] were interchanged.
   
   template<typename real_t>
   int_t sytf2(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * IPIV, bool * BSDV)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      const real_t alpha = (1.0 + std::sqrt(17.0))/8.0;
      const real_t one(1.0);
      const complex<real_t> onec(1.0);
      const real_t zero(0.0);
      complex<real_t> * Ak,* Akm1, d11, d22;
      real_t absakk, colmax;
      int_t info = 0, imax = 0, kp;
      
      for (int_t i = 0; i < n; ++i)
      {
         BSDV[i] = 0;
      }
      
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kstep;
         while (k >= 0)
         {
            kstep = 1;
            Ak = A + ldA*k;
            Akm1 = Ak-ldA;
            
            absakk = std::abs(Ak[k]);
            
            if (k > 0)
            {
               imax = latl::imax(k, Ak, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == zero || std::isnan(absakk))
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
                  int_t jmax;
                  real_t rowmax;
                  
                  jmax = imax+1+ latl::imax(k-imax, A+ldA*(imax+1)+imax, ldA);
                  complex<real_t> * Ajmax = A + ldA*jmax;
                  complex<real_t> * Aimax = A + ldA*imax;
                  rowmax = std::abs(Ajmax[imax]);
                  if (imax > 0)
                  {
                     jmax = latl::imax(imax, Aimax, 1);
                     rowmax = std::max(rowmax, std::abs(Aimax[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (std::abs(Aimax[imax]) >= alpha*rowmax)
                     kp = imax;
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k-kstep+1;
               complex<real_t> * Akinc = A + ldA*kinc;
               complex<real_t> * Akp = A + ldA*kp;
               complex<real_t> r1;
               
               if (kp != kinc)
               {
                  latl::swap(kp, Akinc, 1, Akp, 1);
                  latl::swap(kinc-kp-1, Akinc+kp+1, 1, A+ldA*(kp+1)+kp, ldA);
                  complex<real_t> temp = Akinc[kinc];
                  Akinc[kinc] = Akp[kp];
                  Akp[kp] = temp;
                  if (kstep == 2)
                  {
                     temp = Ak[k-1];
                     Ak[k-1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               
               if (kstep == 1)
               {
                  r1 = onec/Ak[k];
                  latl::syr(uplo, k, -r1, Ak, 1, A, ldA);
                  latl::scal(k, r1, Ak, 1);
               }
               else
               {
                  if (k > 1)
                  {
                     complex<real_t> d12;
                     d12 = Ak[k-1];
                     d22 = Akm1[k-1]/d12;
                     d11 = Ak[k]/d12;
                     d12 = (onec/(d11*d22-onec))/d12;
                     
                     complex<real_t> wk, wkm1;
                     for (int_t j = k-2; j >= 0; --j)
                     {
                        complex<real_t> * Aj = A+ldA*j;
                        wkm1 = d12*(d11*Akm1[j]-Ak[j]);
                        wk = d12*(d22*Ak[j]-Akm1[j]);
                        for (int_t i = j; i >= 0; --i)
                        {
                           Aj[i] -= (Ak[i]*wk+Akm1[i]*wkm1);
                        }
                        Ak[j] = wk;
                        Akm1[j] = wkm1;
                     }
                  }
               }
            }
            if (kstep == 1)
            {
               IPIV[k] = kp;
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
      }
      else
      {
         int_t k = 0, kstep = 1;
         while (k < n)
         {
            kstep = 1;
            Ak = A + ldA*k;
            
            absakk = std::abs(Ak[k]);
            
            if (k < n-1)
            {
               imax = k +1 + latl::imax(n-k-1, Ak+k+1, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == zero || std::isnan(absakk))
            {
               if (info == 0)
                  info = k;
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
                  kp = k;
               else
               {
                  int_t jmax = k + latl::imax(imax-k, Ak+imax, ldA);
                  complex<real_t> * Ajmax = A + ldA*jmax;
                  complex<real_t> *Aimax = A +ldA*imax;
                  real_t rowmax = std::abs(Ajmax[imax]);
                  if (imax < n-1)
                  {
                     jmax = imax +latl::imax(n-imax-1, Aimax+imax+1, 1);
                     rowmax = std::max(rowmax, std::abs(Aimax[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (std::abs(Aimax[imax]) >= alpha*rowmax)
                  {
                     kp = imax;
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k + kstep -1;
               complex<real_t> * Akinc = A+ldA*kinc;
               complex<real_t> * Akp = A+ldA*kp;
               if (kp != kinc)
               {
                  if (kp < n-1)
                     latl::swap(n-kp-1, Akinc+kp+1, 1, Akp+kp+1, 1);
                  latl::swap(kp-kinc-1, Akinc+kinc+1, 1, Akinc+ldA+kp, ldA);
                  complex<real_t> temp = Akinc[kinc];
                  Akinc[kinc] = Akp[kp];
                  Akp[kp] = temp;
                  if (kstep == 2)
                  {
                     temp = Ak[k+1];
                     Ak[k+1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               if (kstep == 1)
               {
                  if (k < n-1)
                  {
                     complex<real_t> r1 = onec/Ak[k];
                     latl::syr(uplo, n-k-1, -r1, Ak+k+1, 1, Ak+ldA+k+1, ldA);
                     latl::scal(n-k-1, r1, Ak+k+1, 1);
                  }
               }
               else
               {
                  if (k < n-2)
                  {
                     complex<real_t> d21, temp, wk, wkp1;
                     complex<real_t> * Akp1 = A+ldA*(k+1);
                     d21 = Ak[k+1];
                     d11 = Akp1[k+1]/d21;
                     d22 = Ak[k]/d21;
                     temp = one/(d11*d22-one);
                     d21 = temp/d21;
                     
                     for (int_t j = k+2; j < n; ++j)
                     {
                        wk = d21*(d11*Ak[j]-Akp1[j]);
                        wkp1 = d21*(d22*Akp1[j]-Ak[j]);
                        complex<real_t> * Aj = A + ldA*j;
                        
                        for (int_t i = j; i< n; ++i)
                        {
                           Aj[i] -= Ak[i]*wk+Akp1[i]*wkp1;
                        }
                        Ak[j] = wk;
                        Akp1[j] = wkp1;
                     }
                  }
               }
            }
            if (kstep == 1)
               IPIV[k] = kp;
            else
            {
               IPIV[k] = kp;
               IPIV[k+1] = kp;
               BSDV[k] = 1;
               BSDV[k+1] = 1;
            }
            
            k += kstep;
         }
      }
      return info;
   }
}

#endif
