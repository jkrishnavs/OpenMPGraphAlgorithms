int64_t T = 0 ;


void outputTriangleCounting(graph *G) {
  printf("\nThe total number of Triangles = %lld\n", T);
}



void triangleCounting(graph *G) {
  inittracking("triangleCounting.csv");
  T = 0;
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
	      edge_t s;
              for(s= G->begin[w]; s < G->begin[w+1]; s++) {
		node_t y = G->node_idx[s];
		if(y == u) {
                  T_private = T_private + 1 ;
		  break;
		}
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
