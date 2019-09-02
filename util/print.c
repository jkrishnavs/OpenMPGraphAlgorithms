#include "print.h"

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
  case CONFIG_FILE_NOT_FOUND:
     printf("Error: Input Config File Not Found. \n");
     break;
  case GRAPH_FILE_NOT_CREATED:
     printf("Error: Unable to create output Graph File. \n");
     break;   
  case OUT_OF_MEMORY:
    printf("Error: Unable to allocate memory.\n");
    break;
  case ARRAY_ACCESS_OUT_OF_BOUNDS:
    printf("Error: Element accessed beyond array Size. \n");
    break;
  case INCORRECT_ARG_LIST:
    printf("Error: Wrong Usage of Arguments: \n Usage: EXECUTABLE");
    int i;
    for(i=0; i< NoOfMsgs;i++) {
      printf("\t%s", msgs[i]);
    }
    printf("\n");
    break;
  case TASKLOOP_NOTENABLED:
    printf("Error: Taskloop not available with the current system\n");
    break;
  case INCORRECT_CONFIG_FILE:
    printf("Error: incorrect Config file, see docs for correct files\n");
    break;
  case UNKNOWN_PROPERTY_ERROR:
    printf("Error: Unknown property in the config file:\n");
    printf("\t%s", msgs[0]);
    printf("\n");
    break;
  default:
    fprintf(stderr,"Unknown error\n");
  }
}
