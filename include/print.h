#ifndef PRINT_DATA_H
#define PRINT_DATA_H
/***
 * This file contains all user interactions,
 * Outputs, messages, error and Warnings.
 **/
#include<stdio.h>
#include "graphEnum.h"

void printTiming(executionSection section, double timeelapsed) {

  switch(section) {
  case GRAPHREAD:
    printf("Reading from File, Time Elapsed =\t%lf ms\n", timeelapsed);
    break;
  case GRAPHWRITE:
    printf("Writing to File, Time Elapsed =\t%lf ms\n", timeelapsed);
    break;
  case REVERSE_EDGE_CREATION:
    printf("Reverse Edge Creation, Time Elapsed =\t%lf ms\n", timeelapsed);
     break;
  case ALGO_KERNEL:
    printf("Algorithm Running Time =\t%lf ms\n", timeelapsed);
    break;
  case OVERALL:
    printf("Overall Running Time =\t%lf ms\n", timeelapsed);
    break;
  default:
    printf("Time Elapsed =\t%lf ms\n",  timeelapsed);
  }

}

void pringMsg(int NoOfMsgs, char**msgs) {
  int i;
  for(i=0; i< NoOfMsgs;i++) {
    printf("\t%s", msgs[i]);
  }
}


void printError(errorCodes code, int NoOfMsgs,const char** msgs) {

  switch(code) {
  case GRAPH_FILE_NOT_FOUND:
     printf("Error: Input Graph File Not Found. \n");
     break;
  case OUT_OF_MEMORY:
    printf("Error: Unable to allocate memory.\n");
    break;
  case ARRAY_ACCESS_OUT_OF_BOUNDS:
    printf("Error: Element accessed beyond array Size. \n");
    break;
  case INCORRECT_ARG_LIST:
    printf("Error: WOrng Usage of Arguments: \n Usage: EXECUTABLE");
    int i;
    for(i=0; i< NoOfMsgs;i++) {
      printf("\t%s", msgs[i]);
    }
    printf("\n");
    break;
  case TASKLOOP_NOTENABLED:
    printf("Error: Taskloop not available with the current system\n");
  default:
    fprintf(stderr,"Unknown error\n");
  }

}

#endif

