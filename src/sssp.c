/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include<unistd.h>
#include "graph.h"
#include "mainFunctions.h"
#include "powerperformacetracking.h"
#include "print.h"
#include "sssp.h"
#include <stdlib.h>

#define NO_OF_ARGS 5



/* 
   Benign non determinism present
   The outer

 */

void sssp(graph *G);

void output(graph* G) {
  outputsssp(G);
}


#define numTimes 7

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  long seed  = 0;
  uint32_t maxLength = 10;
  root = 0;

  int flag = 0;


  if(argc > 3) {
    root  = atoi(argv[3]);
    if(root < 0) flag = 1;
  }

  if(argc > 4) {
    maxLength = atoi(argv[4]);
    if(maxLength < 1) flag = 2;
  }

  if(argc > 5) {
    seed = atoi(argv[5]);
    if(seed < 0) flag = 3;
  }
  

  
  if(flag != 0) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "[root = 0]", "maxLength=10]","[seed = 0]"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }

  
  
  graph* G = readGraph(argv[1], argv[2]);
  // comment for final run.
  assert(G->weights != NULL);
  if(G->weights == NULL) {
    G->weights = (int*) malloc (G->numEdges * sizeof(int));
    /*
      Random generation of edge lengths
    */
    srand(seed);
    
    for(edge_t e = 0; e < G->numEdges; e++) {
      G->weights[e] = (rand()%maxLength + 1);
    }
  }
  dist = (uint32_t*) malloc (G->numNodes * sizeof (uint32_t));
  assert(dist != NULL);
  int i;
  for(i=0;i< numTimes;i++) {
    runKernel(G);
    output(G);
    
    char sssp[50];
    sprintf(sssp, "sssp.%d.csv", i);
    rename("sssp.csv",sssp);

    sleep(5);
  }
  free(dist);
  return 0;
}






inline void kernel(graph *G) {
  sssp(G);
}















