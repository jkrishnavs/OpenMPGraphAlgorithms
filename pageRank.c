/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "include/graph.h"
#include "include/mainFunctions.h"
#include "include/parsegraph.h"
#include "include/print.h"

#define NO_OF_ARGS 4

double e;
double d;
int32_t maxIters;
float* pg_rank;


void pageRank(graph* G) {
  inittracking();
  float eprime = (float) e;
  float diff = 0.0 ;
  int32_t cnt = 0 ;
  double N = 0.0 ;
  
  double* pg_rank_nxt = gm_rt_allocate_double(G.num_nodes(),gm_rt_thread_id());

  cnt = 0 ;
  N = G->numNodes;
  
#pragma omp parallel for
  for (node_t t0 = 0; t0 < G->numNodes; t0 ++) 
    pg_rank[t0] = 1 / N ;

  do
    {
      diff = ((float)(0.000000)) ;
#pragma omp parallel
      {
	float diff_prv = 0.0 ;
	
	diff_prv = ((float)(0.000000)) ;

	
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
	  for (node_t t = 0; t < G->numNodes; t ++) 
            {
	      double val = 0.0 ;
	      double __S1 = 0.0 ;
	      
	      __S1 = ((float)(0.000000)) ;
	      for (edge_t w_idx = G->r_begin[t];w_idx < G->r_begin[t+1] ; w_idx ++) 
                {
		  node_t w = G->r_node_idx [w_idx];
		  __S1 = __S1 + pg_rank[w] / ((double)((G->begin[w+1] - G.begin[w]))) ;
                }
	      val = (1 - d) / N + d * __S1 ;
	      diff_prv = diff_prv +  abs((val - pg_rank[t]))  ;
	      pg_rank_nxt[t] = val ;
            }
#pragma omp atomic
	diff += diff_prv;
      }
      
#pragma omp parallel for
      for (node_t i3 = 0; i3 < G->numNodes; i3 ++) 
	pg_rank[i3] = pg_rank_nxt[i3] ;

        cnt = cnt + 1 ;
    }
    while ((diff > eprime) && (cnt < max));


    gm_rt_cleanup();
    energymonitor__stopprofiling();
}


void output(graph *G) {
  for (int i = 0; i < 4; i++) {
    printf("rank[%d] = %0.9lf\n", i, pg_rank[i]);
  }
  free(pg_rank);
  return true;
}


/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  int flag = 0;
  /* Default values  */
  max_iter = 100;
  e = 0.001;
  d = 0.85;
  
  if(argc > 1) {
    maxIters = atoi(argv[1]);
    if(maxIters <=0) flag = 1;
  }
  
  if(argc > 2) {
    e = atof(argv[2]);
    if(e <= 0) flag = 2;
  }

  if(argc > 3) {
    d = atof(argv[3]);
    if(d <= 0) flag = 3;
  }
  graph* G = parseGraph(argv[1]);
  if(G == NULL)
    flag = 4;
  
  if(flag > 0) {
    const char argList[NO_OF_ARGS] = {" <inputfile>", "[max_iteration=100]", "[eplision=0.001]", "[delta=0.85]" };
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  pg_rank = (float*) malloc (G->numNodes * sizeof(float));
  assert(pg_rank != NULL);
}




















