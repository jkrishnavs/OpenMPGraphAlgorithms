#ifndef POWER_PERFORMANCE_TRACKING_H
#define POWER_PERFORMANCE_TRACKING_H
/******************
 * Add tracking for power and performance here
 * Code will change for different hardwares.
 *
 ******************/
#ifdef BIGLITTE_H
#include "energylib.h"
#endif
void inittracking(char* profName) {
#ifdef BIGLITTE_H
  
  energymonitor__setfilename(profName);
  energymonitor__init(ONLINECORES,0.2);
  energymonitor__trackpoweronly();
  energymonitor__startprofiling();
#endif
}


void pausetracking() {


}


void endtracking() {
#ifdef BIGLITTLE_H
  energymonitor__stopprofiling();
#endif
}


#endif
