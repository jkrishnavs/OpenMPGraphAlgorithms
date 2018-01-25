float __conduct = 0;

int32_t *G_member;

void conduct(graph *G) {
  inittracking();
  for (int i = 0; i < 4; i++) {  
    float m = 0.0 ;
    int32_t __S2 = 0 ;
    int32_t __S3 = 0 ;
    int32_t __S4 = 0 ;
    
#pragma omp parallel
    {
      int32_t __S2_prv = 0 ;
      __S2_prv = 0 ;


#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp for schedule(static)
#endif
      for (node_t u = 0; u < G->numNodes; u ++) 
	if ((G_member[u] == i))
	  {
	    __S2_prv = __S2_prv + (G->begin[u+1] - G->begin[u]) ;
	  }
      
#pragma omp atomic
      __S2 += __S2_prv;
    }

#pragma omp parallel
    {
      int32_t __S3_prv = 0 ;
      __S3_prv = 0 ;
#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp for schedule(static)
#endif
      for (node_t u0 = 0; u0 < G->numNodes; u0 ++) 
	if ((G_member[u0] != i))
	  {
	    __S3_prv = __S3_prv + (G->begin[u0+1] - G->begin[u0]) ;
	  }

#pragma omp atomic
      __S3 += __S3_prv;
    }
    

    
#pragma omp parallel
    {
      int32_t __S4_prv = 0 ;
      __S4_prv = 0 ;

#if defined(PARFOR_GUIDED)   
#pragma omp for schedule(guided, PAR_CHUNKSIZE)
#elif defined(PARFOR_DYNAMIC)
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
#elif defined(TASKLOOP_DEFINED)
#pragma omp taskloop num_tasks(NUM_TASKS)
#else
#pragma omp for schedule(static)
#endif
      for (node_t u1 = 0; u1 < G->numNodes; u1 ++) 
	if ((G_member[u1] == i))
	  {
	    for (edge_t j_idx = G->begin[u1];j_idx < G->begin[u1+1] ; j_idx ++) 
	      {
		node_t j = G->node_idx [j_idx];
	      if ((G_member[j] != i))
		{
		  __S4_prv = __S4_prv + 1 ;
		}
	      }
	  }
#pragma omp atomic
      __S4 += __S4_prv;
    }
    
    m = (float)((__S2 < __S3)?__S2:__S3) ;
  
  }
  endtracking();
  
}
