/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "nodeIntMap.h"
#include "communities.h"

#define NO_OF_ARGS 2


void communities(graph* G);




void output(graph *G) {
  outputCommunities(G);
}


/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  int flag = 0;
  maxItrs = 10;

  if(argc > 2) {
    maxItrs = atoi(argv[2]);
    if(maxItrs <=0) flag = 1;
  }

  graph* G = readGraph(argv[1]);

  if(G == NULL) {
    flag = 2;
  }
  
  if(flag != 0 ) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "[max_iterations = 10]"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  
  initCommunities();
  runKernel(G);
  output(G);
  return 0;
}










inline void kernel(graph *G) {
  communities(G);
}










