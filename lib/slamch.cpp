//
//  slamch.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include <limits>
#include <cctype>

using std::numeric_limits;
using std::toupper;

float slamch_(char &OPT)
{
   const float one=1.0;
   const float zero=0.0;
   const float prec=numeric_limits<float>::epsilon();
   const float rnd=numeric_limits<float>::round_error();
   const float eps=prec*rnd;
   const float maxval=numeric_limits<float>::max();
   const float minval=numeric_limits<float>::min();
   const float small=one/maxval;
   const float sfmin=(small>minval) ? small*(one+eps) : minval;
   const float base=numeric_limits<float>::radix;
   const float num=numeric_limits<float>::digits;
   const float round=numeric_limits<float>::round_style;
   const float min_exp=numeric_limits<float>::min_exponent;
   const float underflow=numeric_limits<float>::min();
   const float max_exp=numeric_limits<float>::max_exponent;
   const float overflow=numeric_limits<float>::max();
   
   const char opt=toupper(OPT);
   
   if(opt=='E')
      return eps;
   else if(opt=='S')
      return sfmin;
   else if(opt=='B')
      return base;
   else if(opt=='P')
      return prec;
   else if(opt=='N')
      return num;
   else if(opt=='R')
      return round;
   else if(opt=='M')
      return min_exp;
   else if(opt=='U')
      return underflow;
   else if(opt=='L')
      return max_exp;
   else if(opt=='O')
      return overflow;
   else
      return zero;
}


