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
#include "pageRank.h"

#define NO_OF_ARGS 5


void output(graph *G) {
  outputPageRank(G);
}


#define numTimes 7

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  int flag = 0;
  /* Default values  */
  maxIters = 100;
  e = 0.01;
  d = 0.85;
  
  if(argc > 3) {
    maxIters = atoi(argv[3]);
    if(maxIters <=0) flag = 1;
  }
  
  if(argc > 4) {
    e = atof(argv[4]);
    if(e <= 0) flag = 2;
  }

  if(argc > 5) {
    d = atof(argv[5]);
    if(d <= 0) flag = 3;
  }


  if(flag != 0) {
    const char* argList[NO_OF_ARGS] = {" <inputfile>", "graphformat.txt","[max_iteration=100]", "[eplision=0.001]", "[delta=0.85]" };
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }

  
  
  graph* G = readGraph(argv[1], argv[2]);
  if(G == NULL)
    flag = 4;
  int i=0;
  createReverseEdge(G);
  for(i=0;i<numTimes;i++) {
    pg_rank = (double*) malloc (G->numNodes * sizeof(double));
    assert(pg_rank != NULL);
    runKernel(G);
    output(G);
    char page[50];
    sprintf(page, "pageRank.%d.csv", i);
    rename("pageRank.csv", page);
    sleep(5);
  }
  //free(pg_rank);
  pg_rank = NULL;
  return 0;
}




inline void kernel(graph* G) {
  pageRank(G);
}












