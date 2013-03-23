//
//  lascl.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 7/2/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lascl_h
#define _lascl_h

/// @file lascl.h Multiplies a matrix by a scalar.

#include <cctype>
#include <cmath>
#include <algorithm>
#include <limits>
#include "latl.h"

namespace LATL
{
   /// @brief  Multiplies real matrix A by the real scalar a/b.
   ///
   /// Multiplication of real matrix A by scalar a/b is done without over/underflow as long as the final
   /// result a*A/b does not over/underflow. The parameter type specifies that
   /// A may be full, upper triangular, lower triangular, upper Hessenberg, or banded.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param type Specifies the type of matrix A.
   ///
   ///        'G':  A is a full matrix.
   ///        'L':  A is a lower triangular matrix.
   ///        'U':  A is an upper triangular matrix.
   ///        'H':  A is an upper Hessenberg matrix.
   ///        'B':  A is a symmetric band matrix with lower bandwidth kl and upper bandwidth ku
   ///              and with the only the lower half stored.
   ///        'Q':  A is a symmetric band matrix with lower bandwidth kl and upper bandwidth ku
   ///              and with the only the upper half stored.
   ///        'Z':  A is a band matrix with lower bandwidth kl and upper bandwidth ku.
   /// @param kl The lower bandwidth of A, used only for banded matrix types B, Q and Z.
   /// @param ku The upper bandwidth of A, used only for banded matrix types B, Q and Z.
   /// @param b The denominator of the scalar a/b.
   /// @param a The numerator of the scalar a/b.
   /// @param m The number of rows of the matrix A. m>=0
   /// @param n The number of columns of the matrix A. n>=0
   /// @param A Pointer to the real matrix A [in/out].
   /// @param ldA The column length of the matrix A.
   
   template<typename real_t>
   int lascl(char type,int_t kl,int_t ku,real_t b,real_t a,int_t m,int_t n,real_t *A,int_t ldA)
   {
      using std::toupper;
      using std::isnan;
      using std::max;
      using std::min;
      using std::numeric_limits;
      const real_t zero=0.0;
      type=toupper(type);
      if((type!='G')&&(type!='L')&&(type!='U')&&(type!='H')&&(type!='B')&&(type!='Q')&&(type!='Z'))
         return -1;
      else if((b==zero)||isnan(b))
         return -4;
      else if(isnan(a))
         return -5;
      else if(m<0)
         return -6;
      else if((n<0)||((type=='B')&&(n!=m))||((type=='Q')&&(n!=m)))
         return -7;
      else if((ldA<m)&&((type=='G')||(type=='L')||(type=='U')||(type=='H')))
         return -9;
      else if(((type=='B')||(type=='Q')||(type=='Z'))&&((kl<0)||(kl>max(m-1,0))))
         return -2;
      else if(((type=='B')||(type=='Q')||(type=='Z'))&&((ku<0)||(ku>max(n-1,0))))
         return -3;
      else if(((type=='B')||(type=='Q'))&&((kl!=ku)))
         return -3;
      else if((type=='B')&&(ldA<kl+1))
         return -9;
      else if((type=='Q')&&(ldA<ku+1))
         return -9;
      else if((type=='Z')&&(ldA<2*kl+ku+1))
         return -9;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t small=numeric_limits<real_t>::min();
      const real_t big=one/small;
      bool done=0;
      real_t c;
      while(!done)
      {
         real_t a1;
         real_t b1=b*small;
         if(b1==b)
         {
            c=a/b;
            done=1;
            a1=a;
         }
         else
         {
            a1=a/big;
            if(a1==a)
            {
               c=a;
               done=1;
               b=one;
            }
            else if((abs(b1)>abs(a))&&(a!=zero))
            {
               c=small;
               done=0;
               b=b1;
            }
            else if(abs(a1)>abs(b))
            {
               c=big;
               done=0;
               a=a1;
            }
            else
            {
               c=a/b;
               done=1;
            }
         }
         if(type=='G')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;i<m;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='L')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=j;i<m;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='U')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<m)&&(i<=j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='H')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<m)&&(i<=j+1);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='B')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<=kl)&&(i<n-j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='Q')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=max(ku-j,0);i<=ku;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='Z')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=max(kl+ku-j,kl);i<=min(2*kl+ku,kl+ku+m-j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
      }
      return 0;
   }
   
   /// @brief  Multiplies complex matrix A by the real scalar a/b.
   ///
   /// Multiplication of complex matrix A by scalar a/b is done without over/underflow as long as the final
   /// result a*A/b does not over/underflow. The parameter type specifies that
   /// A may be full, upper triangular, lower triangular, upper Hessenberg, or banded.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param type Specifies the type of matrix A.
   ///
   ///        'G':  A is a full matrix.
   ///        'L':  A is a lower triangular matrix.
   ///        'U':  A is an upper triangular matrix.
   ///        'H':  A is an upper Hessenberg matrix.
   ///        'B':  A is a symmetric band matrix with lower bandwidth kl and upper bandwidth ku
   ///              and with the only the lower half stored.
   ///        'Q':  A is a symmetric band matrix with lower bandwidth kl and upper bandwidth ku
   ///              and with the only the upper half stored.
   ///        'Z':  A is a band matrix with lower bandwidth kl and upper bandwidth ku.
   /// @param kl The lower bandwidth of A, used only for banded matrix types B, Q and Z.
   /// @param ku The upper bandwidth of A, used only for banded matrix types B, Q and Z.
   /// @param b The denominator of the scalar a/b.
   /// @param a The numerator of the scalar a/b.
   /// @param m The number of rows of the matrix A. m>=0
   /// @param n The number of columns of the matrix A. n>=0
   /// @param A Pointer to the complex matrix A [in/out].
   /// @param ldA The column length of the matrix A.
   
   template<typename real_t>
   int_t lascl(char type,int_t kl,int_t ku,real_t b,real_t a,int_t m,int_t n,complex<real_t> *A,int_t ldA)
   {
      using std::toupper;
      using std::isnan;
      using std::max;
      using std::min;
      using std::numeric_limits;
      const real_t zero=0.0;
      type=toupper(type);
      if((type!='G')&&(type!='L')&&(type!='U')&&(type!='H')&&(type!='B')&&(type!='Q')&&(type!='Z'))
         return -1;
      else if((b==zero)||isnan(b))
         return -4;
      else if(isnan(a))
         return -5;
      else if(m<0)
         return -6;
      else if((n<0)||((type=='B')&&(n!=m))||((type=='Q')&&(n!=m)))
         return -7;
      else if((ldA<m)&&((type=='G')||(type=='L')||(type=='U')||(type=='H')))
         return -9;
      else if(((type=='B')||(type=='Q')||(type=='Z'))&&((kl<0)||(kl>max(m-1,0))))
         return -2;
      else if(((type=='B')||(type=='Q')||(type=='Z'))&&((ku<0)||(ku>max(n-1,0))))
         return -3;
      else if(((type=='B')||(type=='Q'))&&((kl!=ku)))
         return -3;
      else if((type=='B')&&(ldA<kl+1))
         return -9;
      else if((type=='Q')&&(ldA<ku+1))
         return -9;
      else if((type=='Z')&&(ldA<2*kl+ku+1))
         return -9;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t small=numeric_limits<real_t>::min();
      const real_t big=one/small;
      bool done=0;
      real_t c;
      while(!done)
      {
         real_t a1;
         real_t b1=b*small;
         if(b1==b)
         {
            c=a/b;
            done=1;
            a1=a;
         }
         else
         {
            a1=a/big;
            if(a1==a)
            {
               c=a;
               done=1;
               b=one;
            }
            else if((abs(b1)>abs(a))&&(a!=zero))
            {
               c=small;
               done=0;
               b=b1;
            }
            else if(abs(a1)>abs(b))
            {
               c=big;
               done=0;
               a=a1;
            }
            else
            {
               c=a/b;
               done=1;
            }
         }
         if(type=='G')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;i<m;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='L')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=j;i<m;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='U')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<m)&&(i<=j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='H')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<m)&&(i<=j+1);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='B')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<=kl)&&(i<n-j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='Q')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=max(ku-j,0);i<=ku;i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
         else if(type=='Z')
         {
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=max(kl+ku-j,kl);i<=min(2*kl+ku,kl+ku+m-j);i++)
                  A[i]*=c;
               A+=ldA;
            }
         }
      }
      return 0;
   }
}

#endif
