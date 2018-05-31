#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "communities.h"
#include "graphprop.h"

typedef struct graphmap {
  node_t newPos;
  //  node_t revPos;
} graphmap;




void readMap(graph *G, graphmap* gm, const char* filename) {
  int r = 1;
  FILE* f;
  f = fopen(filename, "r");
  node_t i,j;
  node_t id = 0;
  while(id < G->numNodes) {
    r = fscanf(f,"%d %d",&i,&j);
    assert(r != EOF);
    assert(i == id);
    gm[id].newPos = j;
    id++;
  }
  fclose(f);
}


void equivalance(graph *G, graph *newG, graphmap *gm) {

  bool hasEdgeWeight = false;
  if(G->weights != NULL) {
    hasEdgeWeight = true;
  }
  
#pragma omp parallel
  {
    node_t x0;
#pragma omp for schedule(dynamic, 1024)
    for (x0 = 0; x0 < G->numNodes; x0 ++) {
      node_t xNew  = gm[x0].newPos;
      for (edge_t y_idx = G->begin[x0];y_idx < G->begin[x0+1] ; y_idx ++) {
	node_t y = G->node_idx [y_idx];
	node_t yNew = gm[y].newPos;
	bool neighbour = false;
	int weight = 0;
	for(edge_t s = newG->begin[xNew]; s < newG->begin[xNew+1]; s++) {
	  node_t dest = newG->node_idx[s];
	  if(dest == yNew) {
	    neighbour = true;
	    if(hasEdgeWeight)
	      weight = newG->weights[s];
	    break;
	  }  
	}
	assert(neighbour == true);
	if(hasEdgeWeight)
	  assert(weight == G->weights[y_idx]);
      }
    }
  }
}


/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  graph* G = readGraph(argv[1], argv[2]);
  graph* newG = readGraph(argv[3], argv[2]);
  graphmap* gm = (graphmap*) malloc (G->numNodes * sizeof(graphmap));
  readMap(G,gm, argv[4]);
  equivalance(G, newG, gm);
}

inline void kernel(graph *G) {

}

