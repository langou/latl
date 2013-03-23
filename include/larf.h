//
//  larf.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 2/13/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _larf_h
#define _larf_h

/// @file larf.h Applies an elementary reflector to a general rectangular matrix.

#include <algorithm>
#include <cctype>
#include "gemv.h"
#include "ger.h"
#include "latl.h"

namespace LATL
{
   /// @brief Applies a real elementary reflector H to a real m-by-n matrix C.
   ///
   /// The elementary reflector H can be applied on either the left or right, with
   ///
   ///        H = I - tau * v * v'
   /// where tau is a real scalar and v is a real vector.
   /// If tau = 0, then H is taken to be the unit matrix.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the elementary reflector H is applied on the left or right.
   ///
   ///              side='L': form  H * C
   ///              side='R': form  C * H
   /// @param m Number of rows of the matrix C.
   /// @param n Number of columns of the matrix C.
   /// @param[in] v Real vector of containing the elementary reflector.
   ///
   ///              If side='R', v is of length n.
   ///              If side='L', v is of length m.
   /// @param incv Increment of the vector v.
   /// @param tau Value of tau in the representation of H.
   /// @param[in,out] C Real m-by-n matrix.  On exit, C is overwritten with
   ///
   ///                H * C if side='L',
   ///             or C * H if side='R'.
   /// @param ldC Column length of matrix C.  ldC >= m.
   /// @param w Workspace vector (optional).  If used, must of of the following length:
   ///
   ///          n if side='L'
   ///          m if side='R'
   /// otherwise, workspace will be allocated and deallocated internally.

   template<typename real_t>
   int larf(char side, int_t m, int_t n, real_t *v, int_t incv, real_t tau, real_t *C, int_t ldC,real_t *w=NULL)
   {
      const real_t one(1.0);
      const real_t zero(0.0);
      using std::toupper;
      side=toupper(side);
      if((side!='L')&&(side!='R'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(incv==0)
         return -5;
      else if(ldC<m)
         return -8;

      bool allocate=w?0:1;
      if(allocate)
         w=new real_t[(side=='L')?n:m];

      if(side=='L')
      {
         gemv<real_t>('T',m,n,one,C,ldC,v,incv,zero,w,1);
         ger<real_t>(m,n,-tau,v,incv,w,1,C,ldC);
      }
      else
      {
         gemv<real_t>('N',m,n,one,C,ldC,v,incv,zero,w,1);
         ger<real_t>(m,n,-tau,w,1,v,incv,C,ldC);
      }
      if(allocate)
         delete [] w;

      return 0;
   }

   /// @brief Applies a complex elementary reflector H to a real m-by-n matrix C.
   ///
   /// The elementary reflector H can be applied on either the left or right, with
   ///
   ///        H = I - tau * v * v'
   /// where tau is a complex scalar and v is a complex vector.
   /// If tau = 0, then H is taken to be the unit matrix.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the elementary reflector H is applied on the left or right.
   ///
   ///              side='L': form  H * C
   ///              side='R': form  C * H
   /// @param m Number of rows of the matrix C.
   /// @param n Number of columns of the matrix C.
   /// @param[in] v Complex vector of containing the elementary reflector.
   ///
   ///              If side='R', v is of length n.
   ///              If side='L', v is of length m.
   /// @param incv Increment of the vector v.
   /// @param tau Value of tau in the representation of H.
   /// @param[in,out] C Complex m-by-n matrix.  On exit, C is overwritten with
   ///
   ///                H * C if side='L',
   ///             or C * H if side='R'.
   /// @param ldC Column length of matrix C.  ldC >= m.
   /// @param w Workspace vector (optional).  If used, must of of the following length:
   ///
   ///          n if side='L'
   ///          m if side='R'
   /// otherwise, workspace will be allocated and deallocated internally.

   template<typename real_t>
   int larf(char side, int_t m, int_t n, complex<real_t> *v, int_t incv, complex<real_t> tau, complex<real_t> *C, int_t ldC, complex<real_t> *w=NULL)
   {
      const real_t one(1.0);
      const real_t zero(0.0);
      using std::toupper;
      side=toupper(side);
      if((side!='L')&&(side!='R'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(incv==0)
         return -5;
      else if(ldC<m)
         return -8;

      bool allocate=w?0:1;
      if(allocate)
         w=new complex<real_t>[(side=='L')?n:m];
         
      if(side=='L')
      {
         gemv<real_t>('C',m,n,one,C,ldC,v,incv,zero,w,1);
         gerc<real_t>(m,n,-tau,v,incv,w,1,C,ldC);
      }
      else
      {
         gemv<real_t>('N',m,n,one,C,ldC,v,incv,zero,w,1);
         gerc<real_t>(m,n,-tau,w,1,v,incv,C,ldC);
      }
      if(allocate)
         delete [] w;
      return 0;
   }
}
#endif
