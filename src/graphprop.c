/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "graphprop.h"

#define NO_OF_ARGS 2



double avgincident;
double scaledT;

void output(graph *G) {
  printf("Avg Adjecency =  %f \n", avgincident);
  printf("\nAvg ClusterCoeff = %f \n",clusterCoeff);
  printf("Avg Edge Distance = %f \n", aed);
  printf("Sparsity = %.7f \n",sparsityMeasure);
  printf("The Scaled percentange Triangles = %f\n\n", scaledT);
}



/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  omp_set_num_threads(2);
  if(argc < NO_OF_ARGS-1) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " };
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  graph* G = readGraph(argv[1], argv[2]);
  runKernel(G);
  output(G);
  return 0;
}


inline void kernel(graph *G) {
  avgincident = ((double) G->numEdges )/ G->numNodes;
  avgClusterCoeff(G);
  avgEdgeDistance(G);
  diameter(G);
  sparsity(G);
  triangle_counting(G);
  scaledT = ((double) T )/ G->numNodes; 
  sccIndex(G); 
}



















