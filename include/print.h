#ifndef PRINT_DATA_H
#define PRINT_DATA_H
/***
 * This file contains all user interactions,
 * Outputs, messages, error and Warnings.
 **/
#include<stdio.h>
#include "graphEnum.h"
void printTiming(executionSection section, double timeelapsed);
void printMsg(int NoOfMsgs,const char**msgs);
void printError(errorCodes code, int NoOfMsgs,const char** msgs);
#endif

