#define _GNU_SOURCE
#include <syscall.h>
#include <sched.h>
#include "graph.h"
#include "mainFunctions.h"
#include "powerperformacetracking.h"
#include "print.h"
#include <stdlib.h>

#include<unistd.h>

#define NO_OF_ARGS 2

//#define REPEAT 25
#define REPEAT 25

long long iters[8];

struct timeval start, end;

// We define all additional paramemter here
void setaffinity() {
  
  /* #pragma omp parallel
  {
    cpu_set_t newcpu;
    int threadid = omp_get_thread_num();
    CPU_ZERO(&newcpu); 
    CPU_SET ( threadid  , &newcpu) ;
    int __t = sched_setaffinity ( syscall ( SYS_gettid ) , sizeof ( newcpu ) , &newcpu ) ;
    assert(__t == 0);
  }
  */
}






void updateatomicadd(graph *G, int id) {

  printf("The update atomic add %d \n", id);
  char title[50];
  sprintf(title, "updateatomic_%d.csv",id);

  gettimeofday(&start, NULL);
  inittracking(title);
  
  int pf = 0;


#pragma omp parallel
    {
      int flag;
#pragma omp for schedule(dynamic, 1024)
      for (node_t u1 = 0; u1 < G->numNodes; u1 ++) {
	if(u1%2 == 0){
	  for (edge_t u_idx = G->begin[u1];u_idx < G->begin[u1+1] ; u_idx ++) {
            node_t u = G->node_idx [u_idx];
	    if(u1%4 == 0) {
	      for (edge_t w_idx = G->begin[u];w_idx < G->begin[u+1] ; w_idx ++) {
		node_t w = G->node_idx [w_idx];
		flag = 0;
		if (u1%8 == 0) {
		  for(edge_t s = G->begin[u1]; s < G->begin[u1+1]; s++) {
		    node_t y = G->node_idx[s];
		    if(y == w) {
		      flag = 1;
		    }
		  }
		}
	      }
	    }
	  }
	}
	if(u1%2 == 0) {
#pragma omp atomic
	  pf++;
	}
      }	  
    }
    
  endtracking();
  gettimeofday(&end, NULL);

  printf("The pf value is %d \n",pf);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  //free(G_member);
}






#define numTimes 7

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {

  int i;
  setaffinity();
  graph* G = readGraph(argv[1], argv[2]);
  for(i = 0;i< numTimes; i++) {
    printf("Run %d \n", i);
    updateatomicadd(G,i);
    sleep(2);
  }
  return 0;
}

inline void kernel(graph *G) {
  
}


