/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"


#define NO_OF_ARGS 1
int64_t T = 0 ;



void triangle_counting(graph *G) {
  inittracking();
#pragma omp parallel
  {
    int64_t T_private  = 0;
    node_t v;
#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp  for schedule(static)
#endif
    for (v = 0; v < G->numNodes; v ++) {
      edge_t u_idx;
      for (u_idx = G->begin[v]; u_idx < G->begin[v+1]; u_idx ++) {
	node_t u = G->node_idx [u_idx];
	if (u > v) {
	  edge_t w_idx;
	  for (w_idx = G->begin[v]; w_idx < G->begin[v+1]; w_idx ++) {
	    node_t w = G->node_idx [w_idx];
	    if (w > u) {
	      if (isNeighbour(G,w,u)) {
		T_private = T_private + 1 ;
	      }
	    }
	  }
	}
      }
    }
#pragma omp atomic
      T += T_private;
  }
  endtracking();
}

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
  triangle_counting(G);
}



















