/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include <unistd.h>
#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "triangleCounting.h"


#define NO_OF_ARGS 2




void output(graph *G) {
  outputTriangleCounting(G);
}


#define numTimes 7



/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  if(argc < NO_OF_ARGS-1) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> ", "graphformat.txt"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  graph* G = readGraph(argv[1], argv[2]);
  int i;
  for(i=0;i<numTimes;i++) {
    runKernel(G);
    output(G);
    char tri[50];
    sprintf(tri, "triangleCounting.%d.csv", i);
    rename("triangleCounting.csv", tri);
    sleep(5);
  }
  return 0;
}


inline void kernel(graph *G) {
  triangleCounting(G);
}



















