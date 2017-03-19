/* 
 * File:   StopWatch.cpp
 * Author: KjellKod
 * From: https://github.com/KjellKod/StopWatch
 * 
 * Created on 2014-02-07 
 */

#include "StopWatch.h"
StopWatch::StopWatch(){
}

/// @return the elapsed clock cycles since start
clock_t StopWatch::ElapsedUs() const {
   return clock()-mStart;
}

/**
 * Resets the start point
 */
void StopWatch::Restart() {
   mStart = clock();
}


