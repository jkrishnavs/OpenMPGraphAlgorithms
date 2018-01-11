#ifndef POWER_PERFORMANCE_TRACKING_H
#define POWER_PERFORMANCE_TRACKING_H
/******************
 * Add tracking for power and performance here
 * Code will change for different hardwares.
 *
 ******************/
#include "energylib.h"

void inittracking() {
  energymonitor__setfilename("prof.csv");
  energymonitor__init(ONLINECORES,0.2);
  energymonitor__trackpoweronly();
  energymonitor__startprofiling();
}


void pausetracking() {


}


void endtracking() {
  energymonitor__stopprofiling();
}


#endif
