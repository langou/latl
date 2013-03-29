//
//  gtsv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/22/13.
//
//

#ifndef _gtsv_h
#define _gtsv_h

#include "latl.h"

namespace LATL
{

   /// @ingroup DRIV
   
   template< typename real_t>
   int_t GTSV(const int_t n, const int_t nrhs, real_t * const DL, real_t * const D, real_t * const DU, real_t * const B, const int_t ldB)
   {
      if (n < 0)
         return -1;
      if ( nrhs < 0)
         return -2;
      if (ldB < n)
         return -7;
      
      using std::abs;
      if (n == 0 || nrhs == 0)
         return 0;
      
      const real_t zero(0.0);
      real_t fact, temp;
      if (nrhs == 1)
      {
         for (int_t i = 0; i < n-2; ++i)
         {
            if (abs(D[i]) >= abs(DL[i]))
            {
               if ( D[i] != zero)
               {
                  fact = DL[i]/D[i];
                  D[i+1] -=fact*DU[i];
                  B[i+1] -=fact*B[i];
               }
               else
               {
                  return i+1;
               }
               DL[i] = zero;
            }
            else
            {
               fact = D[i]/DL[i];
               D[i] = DL[i];
               temp = D[i+1];
               D[i+1] = DU[i] - fact*temp;
               DL[i] = DU[i+1];
               DU[i+1] = -fact*DL[i];
               DU[i] = temp;
               temp = B[i];
               B[i] = B[i+1];
               B[i+1] = temp - fact*B[i+1];
            }
         }
         if ( n > 1)
         {
            int_t nm2 = n-2;
            if (abs(D[nm2]) >= abs(DL[nm2]))
            {
               if (D[nm2] != zero)
               {
                  fact = DL[nm2]/D[nm2];
                  D[n-1] -=fact*DU[nm2];
                  B[n-1] -=fact*B[nm2];
               }
               else
               {
                  return n;
               }
            }
            else
            {
               fact = D[nm2]/DL[nm2];
               D[nm2] = DL[nm2];
               temp = D[n-1];
               D[n-1] = DU[nm2]-fact*temp;
               DU[nm2] = temp;
               temp = B[nm2];
               B[nm2] = B[n-1];
               B[n-1] = temp - fact*B[n-1];
            }
         }
      }
      else
      {
         real_t * Bj;
         for (int_t i = 0; i < n-2; ++i)
         {
            if (abs(D[i]) >= abs(DL[i]))
            {
               if (D[i] != zero)
               {
                  fact = DL[i]/D[i];
                  D[i+1] -= fact*DU[i];
                  Bj = B;
                  for (int_t j = 0; j < nrhs; ++j)
                  {
                     Bj[i+1] -= fact*Bj[i];
                     Bj += ldB;
                  }
               }
               else
               {
                  return i+1;
               }
               DL[i] = zero;
            }
            else
            {
               fact = D[i]/DL[i];
               D[i] = DL[i];
               temp = D[i+1];
               D[i+1] = DU[i] - fact*temp;
               DL[i] = DU[i+1];
               DU[i+1] = -fact*DL[i];
               DU[i] = temp;
               Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  temp = Bj[i];
                  Bj[i] = Bj[i+1];
                  Bj[i+1] = temp - fact*Bj[i+1];
                  Bj += ldB;
               }
            }
         }
         if (n > 1)
         {
            int_t nm2 = n-2;
            if (abs(D[nm2]) >= abs(DL[nm2]))
            {
               if (D[nm2] != zero)
               {
                  fact = DL[nm2]/D[nm2];
                  D[n-1] -= fact*DU[nm2];
                  Bj = B;
                  for (int_t j = 0; j < nrhs; ++j)
                  {
                     Bj[n-1] -= fact*Bj[nm2];
                     Bj += ldB;
                  }
               }
               else
               {
                  return n-1;
               }
            }
            else
            {
               fact = D[nm2]/DL[nm2];
               D[nm2] = DL[nm2];
               temp = D[n-1];
               D[n-1] = DU[nm2] - fact*temp;
               DU[nm2] = temp;
               Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  temp = Bj[nm2];
                  Bj[nm2] = Bj[n-1];
                  Bj[n-1] = temp - fact*Bj[n-1];
                  Bj += ldB;
               }
            }
         }
         if (D[n-1] == zero)
            return n;
      }
      if (nrhs <= 2)
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            Bj[n-1] = Bj[n-1]/D[n-1];
            if (n > 1)
               Bj[n-2] = (Bj[n-2]-DU[n-2]*Bj[n-1])/D[n-2];
            for (int_t i = n-3; i >= 0; --i)
            {
               Bj[i] = (Bj[i]-DU[i]*Bj[i+1]-DL[i]*Bj[i+2])/D[i];
            }
            Bj += ldB;
         }
      }
      else
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            Bj[n-1] = Bj[n-1]/D[n-1];
            if (n > 1)
               Bj[n-2] = (Bj[n-2] - DU[n-2]*Bj[n-1])/D[n-2];
            for (int_t i = n-3; i >= 0; --i)
            {
               Bj[i] = (Bj[i]-DU[i]*Bj[i+1]-DL[i]*Bj[i+2])/D[i];
            }
            
            Bj += ldB;
         }
      }
   }
}

#endif
