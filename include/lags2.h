//
//  lags2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/27/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lags2_h
#define _lags2_h

/// @file lags2.h Computes 2-by-2 matrices satisfying certain properties.

#include <cmath>
#include "lartg.h"
#include "lasv2.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes 2-by-2 orthogonal matrices satisfying certain properties.
   ///
   /// Matrices U, V and Q, are determined such
   /// that if upper is true then
   ///
   ///        U'*A*Q = U'*( a1 a2 )*Q = ( x  0  )
   ///                    ( 0  a3 )     ( x  x  )
   /// 
   ///        V'*B*Q = V'*( b1 b2 )*Q = ( x  0  )
   ///                    ( 0  b3 )     ( x  x  )
   ///
   /// or if upper is false
   ///
   ///        U'*A*Q = U'*( a1 0  )*Q = ( x  x  )
   ///                    ( a2 a3 )     ( 0  x  )
   ///
   ///        V'*B*Q = V'*( b1 0  )*Q = ( x  x  )
   ///                    ( b2 b3 )     ( 0  x  )
   ///
   /// The rows of the transformed A and B are parallel, where
   ///
   ///        U = (  csu  snu ), V = (  csv snv ), Q = (  csq  snq )
   ///            ( -snu  csu )      ( -snv csv )      ( -snq  csq )
   /// @tparam real_t Floating point type.
   /// @param upper If true, the input matrices are upper triangular; if false, they are lower triangular.
   /// @param a1 First diagonal element of input matrix A.
   /// @param a2 Off-diagaonal element of input matrix A.
   /// @param a3 Second diagonal element of input matrix A.
   /// @param b1 First diagonal element of input matrix B.
   /// @param b2 Off-diagonal element of input matrix B.
   /// @param b3 Second diagonal element of input matrix B.
   /// @param csu Diagonal element of output matrix U.
   /// @param snu Off-diagonal element of output matrix U.
   /// @param csv Diagonal element of output matrix V.
   /// @param snv Off-diagonal element of output matrix V.
   /// @param csq Diagonal element of output matrix Q.
   /// @param snq Off-diagonal element of output matrix Q.
   /// @ingroup SCAL
   
   template<typename real_t>
   void LAGS2(bool upper,real_t a1,real_t a2,real_t a3,real_t b1,real_t b2,real_t b3,real_t &csu,real_t &snu,real_t &csv,real_t &snv,real_t &csq,real_t &snq)
   {
      using std::abs;
      using LATL::LARTG;
      using LATL::LASV2;

      const real_t zero=0.0;
      real_t s1,s2,snr,csr,snl,csl;
      if(upper)
      {
         real_t a=a1*b3;
         real_t d=a3*b1;
         real_t b=a2*b1-a1*b2;
         LASV2(a,b,d,s1,s2,snr,csr,snl,csl);
         if((abs(csl)>=abs(snl))||(abs(csr)>=abs(snr)))
         {
            real_t ua11r=csl*a1;
            real_t ua12=csl*a2+snl*a3;
            real_t vb11r=csr*b1;
            real_t vb12=csr*b2+snr*b3;
            real_t aua12=abs(csl)*abs(a2)+abs(snl)*abs(a3);
            real_t avb12=abs(csr)*abs(b2)+abs(snr)*abs(b3);
            real_t r;
            if((abs(ua11r)+abs(ua12))!=zero)
            {
               if(aua12/(abs(ua11r)+abs(ua12))<=avb12/(abs(vb11r)+abs(vb12)))
                  LARTG(-ua11r,ua12,csq,snq,r);
               else
                  LARTG(-vb11r,vb12,csq,snq,r);
            }
            else
            {
               LARTG(-vb11r,vb12,csq,snq,r);
            }
            csu=csl;
            snu=-snl;
            csv=csr;
            snv=-snr;
         }
         else
         {
            real_t ua21=-snl*a1;
            real_t ua22=-snl*a2+csl*a3;
            real_t vb21=-snr*b1;
            real_t vb22=-snr*b2+csr*b3;
            real_t aua22=abs(snl)*abs(a2)+abs(csl)*abs(a3);
            real_t avb22=abs(snr)*abs(b2)+abs(csr)*abs(b3);
            real_t r;
            if((abs(ua21)+abs(ua22))!=zero)
            {
               if(aua22/(abs(ua21)+abs(ua22))<=avb22/(abs(vb21)+abs(vb22)))
                  LARTG(-ua21,ua22,csq,snq,r);
               else
                  LARTG(-vb21,vb22,csq,snq,r);
            }
            else
            {
               LARTG(-vb21,vb22,csq,snq,r);
            }
            csu=snl;
            snu=csl;
            csv=snr;
            snv=csr;
         }
      }
      else
      {
         real_t a=a1*b3;
         real_t d=a3*b1;
         real_t c=a2*b3-a3*b2;
         LASV2(a,c,d,s1,s2,snr,csr,snl,csl);
         if((abs(csr)>=abs(snr))||(abs(csl)>=abs(snl)))
         {
            real_t ua21=-snr*a1+csr*a2;
            real_t ua22r=csr*a3;
            real_t vb21=-snl*b1+csl+b2;
            real_t vb22r=csl*b3;
            real_t aua21=abs(snr)*abs(a1)+abs(csr)*abs(a2);
            real_t avb21=abs(snl)*abs(b1)+abs(csl)*abs(b2);
            real_t r;
            if((abs(ua21)+abs(ua22r))!=zero)
            {
               if(aua21/(abs(ua21)+abs(ua22r))<=avb21/(abs(vb21)+abs(vb22r)))
                  LARTG(ua22r,ua21,csq,snq,r);
               else
                  LARTG(vb22r,vb21,csq,snq,r);
            }
            else
            {
               LARTG(vb22r,vb21,csq,snq,r);
            }
            csu=csr;
            snu=-snr;
            csv=csl;
            snv=-snl;
         }
         else
         {
            real_t ua11=csr*a1+snr*a2;
            real_t ua12=snr*a3;
            real_t vb11=csl*b1+snl*b2;
            real_t vb12=snl*b3;
            real_t aua11=abs(csr)*abs(a1)+abs(snr)*abs(a2);
            real_t avb11=abs(csl)*abs(b1)+abs(snl)*abs(b2);
            real_t r;
            if((abs(ua11)+abs(ua12))!=zero)
            {
               if(aua11/(abs(ua11)+abs(ua12))<=avb11/(abs(vb11)+abs(vb12)))
                  LARTG(ua12,ua11,csq,snq,r);
               else
                  LARTG(vb12,vb11,csq,snq,r);
            }
            else
            {
               LARTG(vb12,vb11,csq,snq,r);
            }
            csu=snr;
            snu=csr;
            csv=snl;
            snv=csl;
         }
      }
   }
   
   /// @brief Computes 2-by-2 unitary matrices satisfying certain properties.
   ///
   /// Matrices U, V and Q, are determined such
   /// that if upper is true then
   ///
   ///        U'*A*Q = U'*( a1 a2 )*Q = ( x  0  )
   ///                    ( 0  a3 )     ( x  x  )
   /// 
   ///        V'*B*Q = V'*( b1 b2 )*Q = ( x  0  )
   ///                    ( 0  b3 )     ( x  x  )
   ///
   /// or if upper is false
   ///
   ///        U'*A*Q = U'*( a1 0  )*Q = ( x  x  )
   ///                    ( a2 a3 )     ( 0  x  )
   ///
   ///        V'*B*Q = V'*( b1 0  )*Q = ( x  x  )
   ///                    ( b2 b3 )     ( 0  x  )
   ///
   /// The rows of the transformed A and B are parallel, where
   ///
   ///        U = (  csu        snu )
   ///            ( -conj(snu)  csu )
   ///
   ///        V = (  csv        snv )
   ///            ( -conj(snv)  csv ) 
   ///
   ///        Q = (  csq        snq )
   ///            ( -conj(snq)  csq )
   /// @tparam real_t Floating point type.
   /// @param upper If true, the input matrices are upper triangular; if false, they are lower triangular.
   /// @param a1 First diagonal element of input matrix A.
   /// @param a2 Off-diagaonal element of input matrix A.
   /// @param a3 Second diagonal element of input matrix A.
   /// @param b1 First diagonal element of input matrix B.
   /// @param b2 Off-diagonal element of input matrix B.
   /// @param b3 Second diagonal element of input matrix B.
   /// @param csu Diagonal element of output matrix U.
   /// @param snu Off-diagonal element of output matrix U.
   /// @param csv Diagonal element of output matrix V.
   /// @param snv Off-diagonal element of output matrix V.
   /// @param csq Diagonal element of output matrix Q.
   /// @param snq Off-diagonal element of output matrix Q.
   /// @ingroup SCAL
   
   template<typename real_t>
   void LAGS2(bool upper,real_t a1,complex<real_t> a2,real_t a3,real_t b1,complex<real_t> b2,real_t b3,real_t &csu,complex<real_t> &snu,real_t &csv,complex<real_t> &snv,real_t &csq,complex<real_t> &snq)
   {
      using std::real;
      using std::imag;
      using std::abs;
      using std::conj;
      using LATL::LARTG;
      using LATL::LASV2;
      
      const real_t zero=0.0;
      const real_t one=1.0;
      real_t s1,s2,snr,csr,snl,csl;
      complex<real_t> r;
      if(upper)
      {
         real_t a=a1*b3;
         real_t d=a3*b1;
         complex<real_t> b=a2*b1-a1*b2;
         real_t fb=abs(b);
         complex<real_t> d1=(fb==zero)?one:b/fb;
         LASV2(a,fb,d,s1,s2,snr,csr,snl,csl);
         if((abs(csl)>=abs(snl))||(abs(csr)>=abs(snr)))
         {
            real_t ua11r=csl*a1;
            complex<real_t> ua12=csl*a2+d1*snl*a3;
            real_t vb11r=csr*b1;
            complex<real_t> vb12=csr*b2+d1*snr*b3;
            real_t aua12=abs(csl)*abs(a2)+abs(snl)*abs(a3);
            real_t avb12=abs(csr)*abs(b2)+abs(snr)*abs(b3);
            if((abs(ua11r)+abs(ua12))==zero)
               LARTG(-complex<real_t>(vb11r),conj(vb12),csq,snq,r);
            else if((abs(vb11r)+abs(vb12))==zero ) 
               LARTG( -complex<real_t>(ua11r),conj(ua12),csq,snq,r);
            else if(aua12/(abs(ua11r)+abs(ua12))<=avb12/(abs(vb11r)+abs(vb12)))
               LARTG( -complex<real_t>(ua11r),conj(ua12),csq,snq,r);
            else
               LARTG( -complex<real_t>(vb11r),conj(vb12),csq,snq,r);
            csu=csl;
            snu=-d1*snl;
            csv=csr;
            snv=-d1*snr;
         }
         else
         {
            complex<real_t> ua21=-conj(d1)*snl*a1;
            complex<real_t> ua22=-conj(d1)*snl*a2+csl*a3;
            complex<real_t> vb21=-conj(d1)*snr*b1;
            complex<real_t> vb22=-conj(d1)*snr*b2+csr*b3;
            real_t aua22=abs(snl)*abs(a2)+abs(csl)*abs(a3);
            real_t avb22=abs(snr)*abs(b2)+abs(csr)*abs(b3);
            if((abs(ua21)+abs(ua22))==zero)
               LARTG(-conj(vb21),conj(vb22),csq,snq,r);
            else if((abs(vb21)+abs(vb22))==zero)
               LARTG(-conj(ua21),conj(ua22),csq,snq,r);
            else if(aua22/(abs(ua21)+abs(ua22))<=avb22/(abs(vb21)+abs(vb22)))
               LARTG(-conj(ua21),conj(ua22),csq,snq,r);
            else
               LARTG(-conj(vb21),conj(vb22),csq,snq,r);
            csu=snl;
            snu=d1*csl;
            csv=snr;
            snv=d1*csr;
         }
      }
      else
      {
         real_t a=a1*b3;
         real_t d=a3*b1;
         complex<real_t> c=a2*b3-a3*b2;
         real_t fc=abs(c);
         complex<real_t> d1=(fc==zero)?one:c/fc;
         LASV2(a,fc,d,s1,s2,snr,csr,snl,csl);
         if((abs(csr)>=abs(snr))||(abs(csl)>=abs(snl)))
         {
            complex<real_t> ua21=-d1*snr*a1+csr*a2;
            real_t ua22r=csr*a3;
            complex<real_t> vb21=-d1*snl*b1+csl*b2;
            real_t vb22r=csl*b3;
            real_t aua21=abs(snr)*abs(a1)+abs(csr)*abs(a2);
            real_t avb21=abs(snl)*abs(b1)+abs(csl)*abs(b2);
            if((abs(ua21)+abs(ua22r))==zero)
               LARTG(complex<real_t>(vb22r ),vb21,csq,snq,r);
            else if((abs(vb21)+abs(vb22r))==zero)
               LARTG(complex<real_t>(ua22r),ua21,csq,snq,r);
            else if(aua21/(abs(ua21)+abs(ua22r))<=avb21/(abs(vb21)+abs(vb22r)))
               LARTG(complex<real_t>(ua22r),ua21,csq,snq,r);
            else
               LARTG(complex<real_t>(vb22r),vb21,csq,snq,r);
            csu=csr;
            snu=-conj(d1)*snr;
            csv=csl;
            snv=-conj(d1)*snl;
         }
         else
         {
            complex<real_t> ua11=csr*a1+conj(d1)*snr*a2;
            complex<real_t> ua12=conj(d1)*snr*a3;
            complex<real_t> vb11=csl*b1+conj(d1)*snl*b2;
            complex<real_t> vb12=conj(d1)*snl*b3;
            real_t aua11=abs(csr)*abs(a1)+abs(snr)*abs(a2);
            real_t avb11=abs(csl)*abs(b1)+abs(snl)*abs(b2);
            if((abs(ua11)+abs(ua12))==zero)
               LARTG(vb12,vb11,csq,snq,r);
            else if((abs(vb11)+abs(vb12))==zero)
               LARTG(ua12,ua11,csq,snq,r);
            else if(aua11/(abs(ua11)+abs(ua12))<=avb11/(abs(vb11)+abs(vb12)))
               LARTG(ua12,ua11,csq,snq,r);
            else
               LARTG(vb12,vb11,csq,snq,r);
            csu=snr;
            snu=conj(d1)*csr;
            csv=snl;
            snv=conj(d1)*csl;
         }
      }
   }
}

#endif
