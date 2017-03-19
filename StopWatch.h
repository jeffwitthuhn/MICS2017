/* 
 * File:   StopWatch.h
 * Author: KjellKod
 * From: https://github.com/KjellKod/StopWatch
 * 
 * Created on 2014-02-07 
 */


#pragma once
#include <time.h>


class StopWatch {
public:

   StopWatch();

   clock_t ElapsedUs() const;

   void Restart();

private:
   clock_t mStart;
};

