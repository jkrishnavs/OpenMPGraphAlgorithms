/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "powerperformacetracking.h"
#include "print.h"
#include "sssp.h"
#include <stdlib.h>
#include <limits.h>

#define NO_OF_ARGS 4


node_t root;
uint32_t* dist;
uint32_t* len;

/* 
   Benign non determinism present
   The outer

 */

void sssp(graph *G);


void output(graph *G) {

  int counter  = 0;
  
  printf("The shortest path to the following nodes from the root Node %d\n", root);
  for (node_t n = 0; n < G->numNodes && counter < 5; n++) {
    if(dist[n] != UINT_MAX && n != root) {
      printf(" %d -> %d : %ud \n", root, n, dist[n]);
      counter ++;
    }
  }  
}



/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  long seed  = 0;
  uint32_t maxLength = 10;
  root = 0;

  int flag = 0;


  if(argc > 2) {
    root  = atoi(argv[2]);
    if(root < 0) flag = 1;
  }

  if(argc > 3) {
    maxLength = atoi(argv[3]);
    if(maxLength < 1) flag = 2;
  }

  if(argc > 4) {
    seed = atoi(argv[4]);
    if(seed < 0) flag = 3;
  }
  

  
  if(flag != 0) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "[root = 0]", "maxLength=10]","[seed = 0]"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }

  
  
  graph* G = readGraph(argv[1]);
  len = (uint32_t*) malloc (G->numEdges * sizeof(uint32_t));
  /*
    Random generation of edge lengths
   */
  srand(seed);

  for(edge_t e = 0; e < G->numEdges; e++) {
    len[e] = (rand()%maxLength + 1);
  }
  dist = (uint32_t*) malloc (G->numNodes * sizeof (uint32_t));
  assert(dist != NULL);
  runKernel(G);
  output(G);
  return 0;
}






inline void kernel(graph *G) {
  sssp(G);
}















