/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "triangleCounting.h"


#define NO_OF_ARGS 1




void output(graph *G) {
  printf("\nThe total number of Triangles = %lld\n", T);
}



/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  if(argc < NO_OF_ARGS-1) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " };
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  graph* G = readGraph(argv[1]);
  runKernel(G);
  output(G);
  return 0;
}


inline void kernel(graph *G) {
  triangleCounting(G);
}



















