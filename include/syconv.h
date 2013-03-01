//
//  syconv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _syconv_h
#define _syconv_h

/// @file syconv.h

#include "latl.h"

namespace latl
{
   template< typename real_t>
   int_t syconv(const char uplo, const char way, const int_t n, real_t * const A, const int_t ldA, int_t * ipiv, bool * bsdv, real_t * Work)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (way != 'C' && way != 'R' && way != 'c' && way != 'r')
         return -2;
      if ( n < 0)
         return -3;
      if (ldA < n)
         return -5;
      
      if (n == 0)
         return 0;
      
      int_t i = 0, ip;
      const real_t zero(0.0);
      if (uplo == 'U' || uplo == 'u')
      {
         if (way == 'C' || way == 'c')
         {
            i = n-1;
            Work[0] = zero;
            real_t * Ai;
            while (i > 0)
            {
               if (bsdv[i] == 1)
               {
                  Ai = A + ldA*i;
                  Work[i] = Ai[i-1];
                  Ai[i-1] = zero;
                  --i;
               }
               Work[i] = zero;
               --i;
            }
            
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i-1];
                        Aj[i-1] = temp;
                     }
                  }
                  --i;
               }
               --i;
            }
         }
         else
         {
            i = 0;
            while (i < n)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  ++i;
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i-1];
                        Aj[i-1] = temp;
                     }
                  }
               }
               ++i;
            }
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 1)
               {
                  real_t * Ai = A + ldA*i;
                  Ai[i-1] = Work[i];
                  --i;
               }
               --i;
            }
         }
      }
      else
      {
         if (way == 'C' || way == 'c')
         {
            real_t * Ai;
            Work[n-1] = zero;
            while (i < n)
            {
               if ( bsdv[i] == 1 && i < n-1)
               {
                  Ai = A + ldA*i;
                  Work[i] = Ai[i+1];
                  Ai[i+1] = zero;
                  ++i;
               }
               Work[i] = zero;
               ++i;
            }
            i = 0;
            while (i < n)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
                  ++i;
               }
               ++i;
            }
         }
         else
         {
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[i];
                        Aj[i] = Aj[ip];
                        Aj[ip] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  --i;
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[i+1];
                        Aj[i+1] = Aj[ip];
                        Aj[ip] = temp;
                     }
                  }
               }
               --i;
            }
            i = 0;
            while (i < n-1)
            {
               if (bsdv[i] == 1)
               {
                  real_t * Ai = A+ldA*i;
                  Ai[i+1] = Work[i];
                  ++i;
               }
               ++i;
            }
         }
      }
      return 0;
   }
}


#endif
