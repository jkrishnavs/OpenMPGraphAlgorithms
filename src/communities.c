/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include<unistd.h>
#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "communities.h"

#define NO_OF_ARGS 3


void communities(graph* G);




void output(graph *G) {
  outputCommunities(G);
}

#define numTimes 7

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  int flag = 0;
  maxItrs = 10;

  if(argc > 3) {
    maxItrs = atoi(argv[3]);
    if(maxItrs <=0) flag = 1;
  }
  //graph* G = NULL;
  graph* G = readGraph(argv[1], argv[2]);

  if(G == NULL) {
    flag = 2;
  }
  
  if(flag != 0 ) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "[max_iterations = 10]"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }  
  int i;
  for(i=0;i< numTimes;i++){
    initCommunities(G);
    runKernel(G);
    output(G);
    char comm[50];
    sprintf(comm, "communities.%d.csv", i);
    rename("communities.csv", comm);
    sleep(5);
  }
  return 0;
}








inline void kernel(graph *G) {
  communities(G);
}










