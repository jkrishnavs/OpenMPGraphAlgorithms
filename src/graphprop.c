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



double aed;

void avgEdgeDistance(graph *G) {
  long * edgeDistances =  (long*) malloc (sizeof(long) * G->numNodes);
  // long * edgeListSize =  (long*) malloc (sizeof(long) * G->numNodes);

  /* printf("The number of nodes is %d \n", G->numNodes); */
  /* printf("The number of edges is %d \n", G->numEdges); */

    node_t v;

  //double contributions;
  #pragma omp parallel
  {
 
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
      //  edgeListSize[v] = G->begin[v+1] - G->begin[v];
      edgeDistances[v] = 0; 
      for (u_idx = G->begin[v]; u_idx < G->begin[v+1]; u_idx ++) {
	node_t u = G->node_idx [u_idx];
	//	printf("The u is %d the v is %d and dist is %d \n", u,v, abs(u-v));
	edgeDistances[v] += abs(u-v);
      }
    }
  }

  aed = 0.0;

  //  node_t v;
  for(v = 0; v < G->numNodes; v ++) {
    aed += ((double)edgeDistances[v])/G->numEdges;
  }
  
  free(edgeDistances);
}

double clusterCoeff;

void avgClusterCoeff(graph *G) {
  double* localClustering = (double*) malloc (sizeof(double) * G->numNodes);

#pragma omp parallel
  {
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
      localClustering[v] = 0;
      
      for (u_idx = G->begin[v]; u_idx < G->begin[v+1]; u_idx ++) {
	node_t u = G->node_idx [u_idx];
	edge_t w_idx;
	// printf("The node is %d %d \n", v,u);
	for (w_idx = u_idx+1; w_idx < G->begin[v+1]; w_idx ++) {
	  node_t w = G->node_idx [w_idx];
	  //printf("The check for neighbour is between %d and %d \n ",u,w);
	  if (isNeighbour(G,w,u)) {
	    localClustering[v] += 1;
	  }
	  if (isNeighbour(G,u,w)) {
	    localClustering[v] += 1;
	  }
	}
      }
      //      printf("The value of local clustering is %f \n", localClustering[v]);
      
      int neighbours = (int) (G->begin[v+1] - G->begin[v]);
      if(neighbours > 1) {
	localClustering[v] = localClustering[v]/(neighbours * (neighbours -1));
      }
      // printf("The value of local clustering is %f \n", localClustering[v]);
    }
  }
  clusterCoeff = 0;

  node_t v;
  for(v =0;v<G->numNodes;v++) {
    clusterCoeff += localClustering[v];
  }

  clusterCoeff = clusterCoeff/G->numNodes;

  free(localClustering);  
}


void diameter(graph *G) {
  //TODO if required
}

double sparsityMeasure;

void sparsity(graph *G) {
  sparsityMeasure =  ((double)G->numEdges) / (G->numNodes * (G->numNodes-1)); 
}

double sccindex;

void sccIndex(graph *G) {

  double* scclist = (double*) malloc (sizeof(double) * G->numNodes);
  int* visited  = (int*) malloc (sizeof(int) * G->numNodes);
  int* stack = (int*) malloc(sizeof(int) * G->numNodes);
  memset(visited, 0, G->numNodes*sizeof(int));
  
  createReverseEdge(G);

  free(stack);
  free(scclist);
  free(visited);
}


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
	//printf("The edge is from %d to %d \n", v, u);
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
  pausetracking();
}

void output(graph *G) {
  printf("\nAvg ClusterCoeff = %f \n",clusterCoeff);
  printf("Avg Edge Distance = %f \n", aed);
  printf("Sparsity = %f \n",sparsityMeasure);
  printf("The total number of Triangles = %ld\n\n", T);
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
  graph* G = readGraph(argv[1]);
  runKernel(G);
  output(G);
  return 0;
}


inline void kernel(graph *G) {
  avgClusterCoeff(G);
  avgEdgeDistance(G);
  diameter(G);
  sparsity(G);
  triangle_counting(G);
  sccIndex(G); 
}



















