//
//  pbtrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/20/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pbtrf_h
#define _pbtrf_h

/// @file pbtrf.h Computes the Cholesky factorization of a symmetric positive definite band matrix A.

#include "pbtf2.h"
#include "gemm.h"
#include "herk.h"
#include "potf2.h"
#include "syrk.h"
#include "trsm.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the Cholesky factorization of a real symmetric positive definite band matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U' * U   if uplo = 'U'
   ///         A = L * L'   if uplo = 'L'
   ///
   /// where U is upper triangular, U' is the transpose of U, and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo = 'U' or the number of sub-diagonals if uplo = 'L'.  kd >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the upper or lower triangle of the symmetric band matrix A, stored in the first kd+1 rows of the array. On exit, the triangular factor U or L from the Cholesky factorization A = U' * U or A = L * L' of the band matrix A in the same storage format as A.
   /// @param ldAB Column length of the array AB.  ldAB >= kd+1;
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t pbtrf(const char uplo, const int_t n, const int_t kd, real_t * const AB, const int_t ldAB, const int_t nb = 64)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (kd < 0 )
         return -3;
      if (ldAB < kd+1)
         return -5;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb > kd)
         return latl::pbtf2(uplo, n, kd, AB, ldAB);
      else
      {
         const real_t one(1.0);
         const real_t zero(0.0);
         const int_t ldWork = nb+1;
         real_t * const Work = new real_t[ldWork*nb];
         real_t * Worktemp = Work, *ABi = AB;
         int_t i2, i3, ib, info = 0, kld = ldAB-1;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < nb; ++j)
            {
               for (int_t i = 0; i < j; ++i)
               {
                  Worktemp[i] = zero;
               }
               Worktemp += ldWork;
            }
            
            for (int_t i = 0; i < n; i += nb)
            {
               ib = std::min(nb, n-i);
               
               info = latl::potf2(uplo, ib, ABi+kd, kld);
               
               if (info != 0)
               {
                  return info+i;
               }
               if ((i+ib) < n)
               {
                  i2 = std::min(kd-ib, n-i-ib);
                  i3 = std::min(ib, n-i-kd);
                  
                  
                  if (i2 > 0)
                  {
                     latl::trsm('L', 'U', 'T', 'N', ib, i2, one, ABi+kd, kld, ABi+ldAB*ib+kd-ib, kld);
                     latl::syrk('U', 'T', i2, ib, -one, ABi+ldAB*ib+kd-ib, kld, one, ABi+ldAB*ib+kd, kld);
                  }
                  
                  if (i3 > 0)
                  {
                     Worktemp = Work;
                     real_t * ABitemp;
                     for (int_t jj = 0; jj < i3; ++jj)
                     {
                        ABitemp = ABi+ldAB*(jj+kd);
                        for (int_t ii = jj; ii < ib; ++ii)
                        {
                           Worktemp[ii] = ABitemp[ii-jj];
                        }
                        Worktemp += ldWork;
                     }
                     latl::trsm('L', 'U', 'T', 'N', ib, i3, one, ABi+kd, kld, Work, ldWork);
                     
                     if (i2 > 0)
                        latl::gemm('T', 'N', i2, i3, ib, -one, ABi+ldAB*ib+kd-ib, kld, Work, ldWork, one, ABi+ldAB*kd+ib, kld);
                     
                     latl::syrk('U', 'T', i3, ib, -one, Work, ldWork, one, ABi+ldAB*kd+kd, kld);
                     
                     Worktemp = Work;
                     for (int_t jj = 0; jj < i3; ++jj)
                     {
                        ABitemp = ABi+ldAB*(jj+kd);
                        for (int_t ii = jj; ii < ib; ++ii)
                        {
                           ABitemp[ii-jj] = Worktemp[ii];
                        }
                        Worktemp += ldWork;
                     }
                  }
               }
               ABi += ldAB*nb;
            }
         }
         else //if lower
         {
            for (int_t j = 0; j < nb; ++j)
            {
               for (int_t i = j+1; i < nb; ++i)
               {
                  Worktemp[i] = zero;
               }
               Worktemp += ldWork;
            }
            
            for (int_t i = 0; i < n; i += nb)
            {
               ib = std::min(nb, n-i);
               
               info = latl::potf2(uplo, ib, ABi, kld);
               if (info != 0)
               {
                  return info+i;
               }
               if (i+ib < n)
               {
                  i2 = std::min(kd-ib, n-i-ib);
                  i3 = std::min(ib, n-i-kd);
                  
                  if (i2 > 0)
                  {
                     latl::trsm('R', 'L', 'T', 'N', i2, ib, one, ABi, kld, ABi+ib, kld);
                     latl::syrk('L', 'N', i2, ib, -one, ABi+ib, kld, one, ABi+ldAB*ib, kld);
                  }
                  
                  if (i3 > 0)
                  {
                     Worktemp = Work;
                     real_t * ABitemp;
                     for (int_t jj = 0; jj < ib; ++jj)
                     {
                        ABitemp = ABi+ldAB*jj;
                        for (int_t ii = 0; ii <= std::min(jj, i3-1); ++ii)
                        {
                           Worktemp[ii] = ABitemp[kd-jj+ii];
                        }
                        Worktemp += ldWork;
                     }
                     
                     latl::trsm('R', 'L', 'T', 'N', i3, ib, one, ABi, kld, Work, ldWork);
                     if (i2 > 0)
                     {
                        latl::gemm('N', 'T', i3, i2, ib, -one, Work, ldWork, ABi+ib, kld, one, ABi+ldAB*ib+kd-ib, kld);
                     }
                     latl::syrk('L', 'N', i3, ib, -one, Work, ldWork, one, ABi+ldAB*kd, kld);
                     
                     Worktemp = Work;
                     for (int_t jj = 0; jj < ib; ++jj)
                     {
                        ABitemp = ABi+ldAB*jj;
                        for (int_t ii = 0; ii <= std::min(jj, i3-1); ++ii)
                        {
                           ABitemp[kd-jj+ii] = Worktemp[ii];
                        }
                        Worktemp += ldWork;
                     }
                  }
               }
               ABi += ldAB*nb;
            }
         }
         delete [] Work;
      }
   }
   
   /// @brief Computes the Cholesky factorization of a complex Hermitian positive definite band matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U' * U   if uplo = 'U'
   ///         A = L * L'   if uplo = 'L'
   ///
   /// where U is upper triangular, U' is the transpose of U, and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo = 'U' or the number of sub-diagonals if uplo = 'L'.  kd >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the upper or lower triangle of the band matrix A, stored in the first kd+1 rows of the array. On exit, the triangular factor U or L from the Cholesky factorization A = U' * U or A = L * L' of the band matrix A in the same storage format as A.
   /// @param ldAB Column length of the array AB.  ldAB >= kd+1;
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t pbtrf(const char uplo, const int_t n, const int_t kd, complex<real_t> * const AB, const int_t ldAB, const int_t nb = 64)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0 )
         return -2;
      if (kd < 0 )
         return -3;
      if (ldAB < kd+1)
         return -5;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb > kd)
         return latl::pbtf2(uplo, n, kd, AB, ldAB);
      else
      {
         const real_t one(1.0);
         const real_t zero(0.0);
         const complex<real_t> onec(1.0);
         const int_t ldWork = nb+1;
         complex<real_t> * const Work = new complex<real_t>[ldWork*nb];
         complex<real_t> * Worktemp = Work, *ABi = AB;
         int_t i2, i3, ib, info = 0, kld = ldAB-1;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < nb; ++j)
            {
               for (int_t i = 0; i < j; ++i)
               {
                  Worktemp[i] = zero;
               }
               Worktemp += ldWork;
            }
            
            for (int_t i = 0; i < n; i += nb)
            {
               ib = std::min(nb, n-i);
               
               info = latl::potf2(uplo, ib, ABi+kd, kld);
               
               if (info != 0)
               {
                  return info+i;
               }
               if (i+ib < n)
               {
                  i2 = std::min(kd-ib, n-i-ib);
                  i3 = std::min(ib, n-i-kd);
                  
                  if (i2 > 0)
                  {
                     latl::trsm('L', 'U', 'C', 'N', ib, i2, onec, ABi+kd, kld, ABi+ldAB*ib+kd-ib, kld);
                     latl::herk('U', 'C', i2, ib, -one, ABi+ldAB*ib+kd-ib, kld, one, ABi+ldAB*ib+kd, kld);
                  }
                  
                  if (i3 > 0)
                  {
                     Worktemp = Work;
                     complex<real_t> * ABitemp;
                     for (int_t jj = 0; jj < i3; ++jj)
                     {
                        ABitemp = ABi+ldAB*(jj+kd);
                        for (int_t ii = jj; ii < ib; ++ii)
                        {
                           Worktemp[ii] = ABitemp[ii-jj];
                        }
                        Worktemp += ldWork;
                     }
                     latl::trsm('L', 'U', 'C', 'N', ib, i3, onec, ABi+kd, kld, Work, ldWork);
                     
                     if (i2 > 0)
                        latl::gemm('C', 'N', i2, i3, ib, -onec, ABi+ldAB*ib+kd-ib, kld, Work, ldWork, onec, ABi+ldAB*kd+ib, kld);
                     
                     latl::herk('U', 'C', i3, ib, -one, Work, ldWork, one, ABi+ldAB*kd+kd, kld);
                     
                     Worktemp = Work;
                     for (int_t jj = 0; jj < i3; ++jj)
                     {
                        ABitemp = ABi+ldAB*(jj+kd);
                        for (int_t ii = jj; ii < ib; ++ii)
                        {
                           ABitemp[ii-jj] = Worktemp[ii];
                        }
                        Worktemp += ldWork;
                     }
                  }
               }
               ABi += ldAB*nb;
            }
         }
         else
         {
            for (int_t j = 0; j < nb; ++j)
            {
               for (int_t i = j+1; i < nb; ++i)
               {
                  Worktemp[i] = zero;
               }
               Worktemp += ldWork;
            }
            
            for (int_t i = 0; i < n; i += nb)
            {
               ib = std::min(nb, n-i);
               
               info = latl::potf2(uplo, ib, ABi, kld);
               if (info != 0)
               {
                  return info+i;
               }
               if (i+ib < n)
               {
                  i2 = std::min(kd-ib, n-i-ib);
                  i3 = std::min(ib, n-i-kd);
                  
                  if (i2 > 0)
                  {
                     latl::trsm('R', 'L', 'C', 'N', i2, ib, onec, ABi, kld, ABi+ib, kld);
                     latl::herk('L', 'N', i2, ib, -one, ABi+ib, kld, one, ABi+ldAB*ib, kld);
                  }
                  
                  if (i3 > 0)
                  {
                     Worktemp = Work;
                     complex<real_t> * ABitemp;
                     for (int_t jj = 0; jj < ib; ++jj)
                     {
                        ABitemp = ABi+ldAB*jj;
                        for (int_t ii = 0; ii <= std::min(jj, i3-1); ++ii)
                        {
                           Worktemp[ii] = ABitemp[kd-jj+ii];
                        }
                        Worktemp += ldWork;
                     }
                     
                     latl::trsm('R', 'L', 'C', 'N', i3, ib, onec, ABi, kld, Work, ldWork);
                     if (i2 > 0)
                     {
                        latl::gemm('N', 'C', i3, i2, ib, -onec, Work, ldWork, ABi+ib, kld, onec, ABi+ldAB*ib+kd-ib, kld);
                     }
                     latl::herk('L', 'N', i3, ib, -one, Work, ldWork, one, ABi+ldAB*kd, kld);
                     
                     Worktemp = Work;
                     for (int_t jj = 0; jj < ib; ++jj)
                     {
                        ABitemp = ABi+ldAB*jj;
                        for (int_t ii = 0; ii <= std::min(jj, i3-1); ++ii)
                        {
                           ABitemp[kd-jj+ii] = Worktemp[ii];
                        }
                        Worktemp += ldWork;
                     }
                  }
               }
               ABi += ldAB*nb;
            }
         }
         delete [] Work;
      }
   }
}


#endif