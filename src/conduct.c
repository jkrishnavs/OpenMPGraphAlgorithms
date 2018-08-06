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
#include "conduct.h"
#include <float.h>

#define NO_OF_ARGS 2



void output(graph *G) {
  outputConduct(G);
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
  G_member = (int32_t*) malloc (G->numNodes * sizeof(int32_t));  
  

  for(i=0;i< numTimes;i++) {
    srand(0);
#pragma parallel for 
    for (int i = 0; i < G->numNodes; i++) 
      G_member[i] = 0;
    
    for (int i = 0; i < G->numNodes; i++) {
      int32_t r = rand() % 100;
      if (r < 10)
	G_member[i] = 0;  // 10%
      else if (r < (10 + 20))
	G_member[i] = 1;  // 20%
      else if (r < (10 + 20 + 30))
	G_member[i] = 2;  // 30%
      else
	G_member[i] = 3;  // 40%
    }  
    runKernel(G);
    output(G);
    char cond[50];
    sprintf(cond, "conduct.%d.csv", i);
    rename("conduct.csv", cond);
    sleep(5);
  }
  free(G_member);
  return 0;
}

inline void kernel(graph *G) {
  conduct(G);
}



















