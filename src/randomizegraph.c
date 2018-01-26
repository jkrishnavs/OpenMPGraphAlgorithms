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

#define NO_OF_ARGS 2

typedef struct graphmap{
  node_t orig;
  node_t random;
  node_t revPos;
}graphmap;

graph* randomizeGraph(graph *G) {
  graphmap* gm = (graphmap*) malloc (G->numNodes * sizeof(graphmap));

  node_t i;
  for(i = 0; i<G->numNodes;i++) {
    gm[i].orig = i;
    gm[i].random = NIL_NODE;
    gm[i].revPos = NIL_NODE;
  }

  srand(time(NULL));
  
  for(i = 0;i< G->numNodes; i++) {
    if(gm[i].random == NIL_NODE) {
      node_t r = (node_t)(rand()% G->numNodes); 
      if(gm[r].random == NIL_NODE) {
	gm[r].random = i;
	gm[i].random = r;

	gm[i].revPos = r;
	gm[r].revPos = i;
      } else {
	// replace new r with old r.
	node_t torep = gm[r].random;
	gm[r].random = i;
	gm[i].random = torep;
	
	gm[i].revPos = r;
	gm[torep].revPos = i; 
      }
    }
  }

  /* Create the copy */
  graph* newG = createGraph();

  newG->numNodes = G->numNodes;
  newG->numEdges = G->numEdges;
  newG->begin = (edge_t*) malloc (sizeof(edge_t) * (newG->numNodes+1));
  assert(newG->begin != NULL);
  newG->node_idx = (node_t*) malloc(sizeof(node_t) * newG->numEdges);

  edge_t edgepos = 0;

  newG->begin[0] = 0;
  for(i=0;i< newG->numNodes; i++) {
    node_t randPos = gm[i].random;
    newG->begin[i+1] = newG->begin[i] + (G->begin[randPos+1] - G->begin[randPos]);
    edge_t st = newG->begin[i];
    edge_t ed = newG->begin[i];
    for(edge_t e = G->begin[randPos]; e < G->begin[randPos+1]; e++) {
      node_t end = G->node_idx[e];
      node_t endRand = gm[end].revPos;
      edge_t e = st;
      for(; e< ed ; e++) {
	if(endRand < newG->node_idx[e])
	  break;
      }
      
      for(edge_t sf = ed; sf >= e; sf --) {
	newG->node_idx[sf+1] = newG->node_idx[sf];  
      }      
      newG->node_idx[e] = endRand;
      ed++;
    }
    assert(ed == newG->begin[i+1]);
  }
  return newG;
}

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  omp_set_num_threads(2);
  if(argc < NO_OF_ARGS-1) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "<outputfile>"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  graph* G = readGraph(argv[1]);
  graph* newG = randomizeGraph(G);
  deleteGraph(G);
  writeGraph(newG, argv[2]);
  return 0;
}



void kernel(graph *G) {


}
















