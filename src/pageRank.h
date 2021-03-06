#include <math.h>
double e;
double d;
int32_t maxIters;
double* pg_rank;



void outputPageRank(graph *G) {
  for (int i = 0; i < 10 && i < G->numNodes; i++) {
    printf("rank[%d] = %0.9lf\n", i, pg_rank[i]);
  }
  free(pg_rank);
}


void pageRank(graph* G) {
  inittracking("pageRank.csv");
  double diff = 0.0 ;
  int32_t cnt = 0 ;
  double N = 0.0 ;
  
  double* pg_rank_nxt = (double*) malloc (G->numNodes * sizeof(double));

  cnt = 0 ;
  N = G->numNodes;
  
#pragma omp parallel for
  for (node_t t0 = 0; t0 < G->numNodes; t0 ++) 
    pg_rank[t0] = 1 / N ;

  do
    {
      diff = 0.0 ;
#pragma omp parallel
      {
	double diff_prv = 0.0 ;
	

	
#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
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
		  __S1 = __S1 + pg_rank[w] / ((double)((G->begin[w+1] - G->begin[w]))) ;
                }
	      val = (1 - d) / N + d * __S1 ;
	      diff_prv = diff_prv +  fabs((val - pg_rank[t]))  ;
	      pg_rank_nxt[t] = val ;
            }
#pragma omp atomic
	diff += diff_prv;
#pragma omp single
	{
	  double * temp = pg_rank_nxt;
	  pg_rank_nxt  = pg_rank;
	  pg_rank = temp;
	}
      }
      cnt = cnt + 1 ;
    } while ((diff > e) && (cnt < maxIters));

  printf("The number if iterations is %d \n",cnt);

  free(pg_rank_nxt);
  endtracking();
}
