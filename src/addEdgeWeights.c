/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include <time.h>
#include <stdlib.h>

#define NO_OF_ARGS 3

typedef struct graphmap{
  node_t orig;
  node_t random;
  node_t revPos;
}graphmap;

int maxLength = 100;
int seed = 0;


graph* addweights(graph *G) {
  if(G->weights != NULL)
    return G;

  G->weights = (int*) malloc (G->numEdges * sizeof(int));
  /*
    Random generation of edge lengths
  */
  srand(seed);
  
  for(edge_t e = 0; e < G->numEdges; e++) {
    G->weights[e] = (rand()%maxLength + 1);
  }
  
  return G;
}



/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  graph* G = readGraph(argv[1], argv[2]);
  graph* newG = addweights(G);
  writeGraph(newG, argv[3]);
  return 0;
}



void kernel(graph *G) {


}
















