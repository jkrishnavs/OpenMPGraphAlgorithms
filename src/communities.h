#include "nodeIntMap.h"

node_t* comm = NULL;

int maxItrs;


void outputCommunities(graph *G) {
  // print output.
  node_t commList[10];
  int commCount[10];
  int i;
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

}


void initCommunities(graph *G) {
  
  comm = (node_t*) malloc (G->numNodes * sizeof(node_t));
  assert(comm != NULL);

}


void communities(graph* G) {
  inittracking("communities.csv");
  bool finished = false ;
  
  finished = true ;
  
#if defined(PARFOR_GUIDED)   
#pragma omp parallel for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp parallel for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp parallel
  {
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp parallel for schedule(static)
#endif
    for (node_t x = 0; x < G->numNodes; x ++) 
      comm[x] = x ;
#if defined(TASKLOOP_DEFINED)
  }
#endif

  int itrs = 0;
  
  do
    {
      finished = true ;


#pragma omp parallel
      {
	nodeIntMap *map;
	map = NULL;
	map = initNodeIntMap(map, 32, 0);
	node_t x0;
#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp  for schedule(static)
#endif
	for (x0 = 0; x0 < G->numNodes; x0 ++) {
	  map = reinitNodeIntMap(map, G->begin[x0+1] - G->begin[x0], 0);
	  for (edge_t y_idx = G->begin[x0];y_idx < G->begin[x0+1] ; y_idx ++) {
	    node_t y = G->node_idx [y_idx];
	    node_t source;
#pragma omp atomic read
	    source = comm[y] ;
	    changeValue(map, source, 1);
	  }
	  node_t maxVal = mapMaxValueKey(map); 
	  if ( comm[x0] != maxVal && maxVal != NIL_NODE) {
#pragma omp atomic write
	    comm[x0] = maxVal;
	    finished = false ;
	  }
	}
	closeNodeIntMap(map);
      }
      itrs++;
    } while ( !finished && maxItrs > itrs);
  
  printf( "Iterations %d \n", itrs);
  endtracking();
}
