//
//  timer.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _timer_h
#define _timer_h

namespace LATL
{
  class Timer
  {
  public:
    Timer();
    ~Timer();
    void Start();
    void Stop();
    double Time();
    void Clear();
  private:
    double StopWatch;
    double Total;
    double Clock();
  };
}
#endif

