/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "include/graph.h"
#include "include/mainFunctions.h"
#include "include/parsegraph.h"
#include "include/print.h"
#include "include/powerandperformancetracking.h"
#include "includenodeIntMap.h"
#define NO_OF_ARGS 1


node_t* comm = NULL;


void communities(graph* G) {
  inittracking();
  bool finished = false ;
  
  finished = true ;
  
#if defined(PARFOR_GUIDED)   
#pragma omp parallel for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp parallel for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP)
#ifndef TASKLOOP_DEFINED
  printError(TASKLOOP_NOTENABLED);
  return
#endif
#pragma omp parallel taskloop
#else
#pragma omp parallel for schedule(static)
#endif
    for (node_t x = 0; x < G.num_nodes(); x ++) 
      comm[x] = x ;
  
  do
    {
      finished = true ;


#pragma omp parallel
      {
	nodeIntMap map;
	initNodeIntMap(map,0);
	
#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP)
#ifndef TASKLOOP_DEFINED
	printError(TASKLOOP_NOTENABLED);
	return
#endif
#pragma omp taskloop
#else
#pragma omp  for schedule(static)
#endif
	  for (node_t x0 = 0; x0 < G.num_nodes(); x0 ++) {
	    reinitNodeIntMap(map, G->begin[x0+1] - G->begin[x0], 0);
	    for (edge_t y_idx = G->begin[x0];y_idx < G->begin[x0+1] ; y_idx ++) {
	      node_t y = G->node_idx [y_idx];
	      node_t source;
	      source = comm[y] ;
	      changeValue(map, source, 1);
	    }
	    node_t maxVal = mapMaxValueKey(map); 
	    if ( comm[x0] != maxVal) {
	      comm[x0] = maxVal;
	      finished = false ;
	    }
	  }
	closeNodeIntMap(map);
      }
    } while ( !finished);
  
  
  endtracking();
}


void output(graph *G) {
  // print output.
  node_t commList[10];
  int commCount[10];
  int i,j;
  for(i=0;i<10; i++) {
    commList[i] = NIL_NODE;
    commCount[i] = 0;
  }
  int found;
  int curIndex = 0;
  int t;
  for (t = 0; t < G->numNodes; t++) {
    found  = -1;
    for(i = 0; i<10;i++) {
      if(comm[t] == commList[i]) {
	found = i;
	break;
      }
    }
    if(found != -1) {
      commCount[found]++;
    } else if(curIndex < 10) {
      commCount[curIndex] = 1;
      commList[curIndex] = comm[t];
      curIndex++;
    }    
  }
  printf("Community\t#Nodes\t\t(Showing max 10 entries)\n");
  for (i=0; i<10; i++) {
    if(commList[i] != NIL_NODE)
      printf("%d\t\t%d\n", commList[i], commCount[i]);
  }
  free(comm);
  comm = NULL;
  return true;
}


/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  if(argc < NO_OF_ARGS-1) {
    const char argList[NO_OF_ARGS] = {" <inputfile> " };
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  graph* G = parseGraph(argv[1]);
  comm = (node_t*) malloc (G->numNodes * sizeof(node_t));
  assert(comm != NULL);
}




















