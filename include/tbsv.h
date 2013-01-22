//
//  tbsv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _tbsv_h
#define _tbsv_h

/// @file tbsv.h Solves banded triangular system of equations.

#include <algorithm>
#include <cctype>
#include "latl.h"

namespace latl
{
   /// @brief Solves banded real triangular system of equations.
   /// 
   /// For a banded real upper or lower triangular matrix A and real vector x,
   ///
   ///       A*x=b  or  A'*x=b
   ///
   /// is solved for x.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose of A is to be used or not:
   ///
   ///        if trans = 'N' or 'n' then x := A*x
   ///        if trans = 'T' or 't' then x := A'*x
   ///        if trans = 'C' or 'c' then x := A'*x
   ///
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param k Specifies the number of super/sub-diagonals of the matrix A.  k >= 0
   /// @param A Pointer to banded triangular n-by-n matrix A.  The bands of A are stored as rows,
   /// while preserving the columns of A.  If uplo = 'U' or 'u' then only the upper triangular part of A 
   /// is referenced and the lower part is not referenced; the diagonal is stored on row k, the first super-diagonal
   /// in row k-1 starting in column 1, the second super-diagonal in row k-2 starting in column 2, and so on.  
   /// As an example, consider the following upper triangular matrix with n=5 and k=2.  On the left is the matrix in
   /// standard upper triangular form, and on the right is the same matrix in upper triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ . i e b . ]        [ . d e f g ]
   ///        [ . . j f c ]        [ h i j k l ]
   ///        [ . . . k g ]        
   ///        [ . . . . l ]
   /// If uplo = 'L' or 'l' then only the lower triangular part of A is referenced and the upper part is not referenced;
   /// the diagonal is stored in row zero, the first sub-diagonal in row 1, and second sub-diagonal in row 2, and so on.  
   /// As an example, consider the following lower triangular matrix with n=5 and k=2.  On the left is the matrix in standard 
   /// lower triangular form, and on the right in the same matrix in lower triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h . . . . ]        [ h i j k l ]
   ///        [ d i . . . ]        [ d e f g . ]
   ///        [ a e j . . ]        [ a b c . . ]
   ///        [ . b f k . ]        
   ///        [ . . c g l ]
   /// @param ldA Column length of matrix A.  ldA>=k+1.
   /// @param x Pointer to real vector x.  On entry, x contains the right hand side vector b.
   /// On exit, x is overwritten with the solution vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup SOLV

   template <typename real_t>
   int tbsv(char uplo, char trans, char diag, int_t n, int_t k, real_t *A, int_t ldA, real_t *x, int_t incx)
   {
      using std::min;
      using std::max;
      using std::toupper;
      
      int_t i,j,ix,jx,kx,i0,in;
      const int_t izero=0;
      bool nounit;

      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if((diag!='U')&&(diag!='N'))
         return -3;
      else if(n<0)
         return -4;
      else if(k<0)
         return -5;
      else if(ldA<k+1)
         return -7;
      else if(incx==0)
         return -9;
      else if(n==0)
         return 0;
      
      nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  if(nounit)
                     x[j]=x[j]/A[k];
                  i0=max(izero,j-k);
                  for(i=j-1;i>=i0;i--)
                     x[i]-=x[j]*A[k-j+i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]=x[j]/A[0];
                  in=min(n-1,j+k);
                  for(i=j+1;i<=in;i++)
                     x[i]-=x[j]*A[i-j];
                  A+=ldA;
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  for(i=i0;i<j;i++)
                     x[j]-=x[i]*A[k-j+i];
                  if(nounit)
                     x[j]=x[j]/A[k];
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  in=min(n-1,j+k);
                  for(i=in;i>j;i--)
                     x[j]-=x[i]*A[i-j];
                  if(nounit)
                     x[j]=x[j]/A[0];
               }
            }
         }
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]=x[jx]/A[k];
                  i0=max(izero,j-k);
                  for(i=j-1;i>=i0;i--)
                  {
                     ix-=incx;
                     x[ix]-=x[jx]*A[k-j+i];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[jx]=x[jx]/A[0];
                  ix=jx;
                  in=min(n-1,j+k);
                  for(i=j+1;i<=in;i++)
                  {
                     ix+=incx;
                     x[ix]-=x[jx]*A[i-j];
                  }
                  jx+=incx;
                  A+=ldA;
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  ix=kx+i0*incx;
                  for(i=i0;i<j;i++)
                  {
                     x[jx]-=x[ix]*A[k-j+i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[k];
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  A-=ldA;
                  in=min(n-1,j+k);
                  ix=kx-(n-in)*incx;
                  for(i=in;i>j;i--)
                  {
                     x[jx]-=x[ix]*A[i-j];
                     ix-=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[0];
               }
            }
         }
      }
      return 0;
   }
   
   /// @brief Solves banded complex triangular system of equations.
   /// 
   /// For a banded complex upper or lower triangular matrix A and real vector x,
   ///
   ///       A*x=b  or  A'*x=b  or A.'*x=b
   ///
   /// is solved for x.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose or conjugate transpose of A is to be used:
   ///
   ///        if trans = 'N' or 'n' then x := A*x
   ///        if trans = 'T' or 't' then x := A.'*x
   ///        if trans = 'C' or 'c' then x := A'*x
   ///
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param k Specifies the number of super/sub-diagonals of the matrix A.  k >= 0
   /// @param A Pointer to banded triangular n-by-n matrix A.  The bands of A are stored as rows,
   /// while preserving the columns of A.  If uplo = 'U' or 'u' then only the upper triangular part of A 
   /// is referenced and the lower part is not referenced; the diagonal is stored on row k, the first super-diagonal
   /// in row k-1 starting in column 1, the second super-diagonal in row k-2 starting in column 2, and so on.  
   /// As an example, consider the following upper triangular matrix with n=5 and k=2.  On the left is the matrix in
   /// standard upper triangular form, and on the right is the same matrix in upper triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ . i e b . ]        [ . d e f g ]
   ///        [ . . j f c ]        [ h i j k l ]
   ///        [ . . . k g ]        
   ///        [ . . . . l ]
   /// If uplo = 'L' or 'l' then only the lower triangular part of A is referenced and the upper part is not referenced;
   /// the diagonal is stored in row zero, the first sub-diagonal in row 1, and second sub-diagonal in row 2, and so on.  
   /// As an example, consider the following lower triangular matrix with n=5 and k=2.  On the left is the matrix in standard 
   /// lower triangular form, and on the right in the same matrix in lower triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h . . . . ]        [ h i j k l ]
   ///        [ d i . . . ]        [ d e f g . ]
   ///        [ a e j . . ]        [ a b c . . ]
   ///        [ . b f k . ]        
   ///        [ . . c g l ]
   /// @param ldA Column length of matrix A.  ldA>=k+1.
   /// @param x Pointer to complex vector x.  On entry, x contains the right hand side vector b.
   /// On exit, x is overwritten with the solution vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup SOLV

   template <typename real_t>
   int tbsv(char uplo, char trans, char diag, int_t n, int_t k, complex<real_t> *A, int_t ldA, complex<real_t> *x, int_t incx)
   {
      using std::conj;
      using std::min;
      using std::max;
      using std::toupper;
      
      int_t i,j,ix,jx,kx,i0,in;
      bool nounit;
      const int_t izero=0;

      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if((diag!='U')&&(diag!='N'))
         return -3;
      else if(n<0)
         return -4;
      else if(k<0)
         return -5;
      else if(ldA<k+1)
         return -7;
      else if(incx==0)
         return -9;
      else if(n==0)
         return 0;
      
      nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  if(nounit)
                     x[j]=x[j]/A[k];
                  i0=max(izero,j-k);
                  for(i=j-1;i>=i0;i--)
                     x[i]-=x[j]*A[k-j+i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]=x[j]/A[0];
                  in=min(n-1,j+k);
                  for(i=j+1;i<=in;i++)
                     x[i]-=x[j]*A[i-j];
                  A+=ldA;
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  for(i=i0;i<j;i++)
                     x[j]-=x[i]*A[k-j+i];
                  if(nounit)
                     x[j]=x[j]/A[k];
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  in=min(n-1,j+k);
                  for(i=in;i>j;i--)
                     x[j]-=x[i]*A[i-j];
                  if(nounit)
                     x[j]=x[j]/A[0];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  for(i=i0;i<j;i++)
                     x[j]-=x[i]*conj(A[k-j+i]);
                  if(nounit)
                     x[j]=x[j]/conj(A[k]);
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  in=min(n-1,j+k);
                  for(i=in;i>j;i--)
                     x[j]-=x[i]*conj(A[i-j]);
                  if(nounit)
                     x[j]=x[j]/conj(A[0]);
               }
            }
         }      
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]=x[jx]/A[k];
                  i0=max(izero,j-k);
                  for(i=j-1;i>=i0;i--)
                  {
                     ix-=incx;
                     x[ix]-=x[jx]*A[k-j+i];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[jx]=x[jx]/A[0];
                  ix=jx;
                  in=min(n-1,j+k);
                  for(i=j+1;i<=in;i++)
                  {
                     ix+=incx;
                     x[ix]-=x[jx]*A[i-j];
                  }
                  jx+=incx;
                  A+=ldA;
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  ix=kx+i0*incx;
                  for(i=i0;i<j;i++)
                  {
                     x[jx]-=x[ix]*A[k-j+i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[k];
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  A-=ldA;
                  in=min(n-1,j+k);
                  ix=kx-(n-in)*incx;
                  for(i=in;i>j;i--)
                  {
                     x[jx]-=x[ix]*A[i-j];
                     ix-=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[0];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  i0=max(izero,j-k);
                  ix=kx+i0*incx;
                  for(i=i0;i<j;i++)
                  {
                     x[jx]-=x[ix]*conj(A[k-j+i]);
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/conj(A[k]);
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  A-=ldA;
                  in=min(n-1,j+k);
                  ix=kx-(n-in)*incx;
                  for(i=in;i>j;i--)
                  {
                     x[jx]-=x[ix]*conj(A[i-j]);
                     ix-=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/conj(A[0]);
               }
            }
         }      
      }
      return 0;
   }
}
#endif
