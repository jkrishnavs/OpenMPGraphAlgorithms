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

#define REPEAT 15

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


void train7(graph *G, int id) {

  printf("The train 7 Array %d \n", id);
  node_t * G_member = (int*)malloc (G->numNodes * sizeof(int));
  srand(0);
  int i;
  for(i = 0; i< G->numNodes; i++) {
    G_member[i] = rand() % G->numNodes;
  }
  char title[50];
  sprintf(title, "train7_%d.csv",id);

  gettimeofday(&start, NULL);
  inittracking(title);
  int tShared = 0;

  for(int abc=0; abc < REPEAT; abc ++) {
    
#pragma omp parallel
    {
      int threadid = omp_get_thread_num();
      int t = 0;
#pragma omp for schedule(dynamic, 1024)
      for (node_t u1 = 0; u1 < G->numNodes; u1 ++)
	for (edge_t j_idx = G->begin[u1];j_idx < G->begin[u1+1] ; j_idx ++) {
	  node_t j = G->node_idx [j_idx];	  
	  for (edge_t k_idx = G->begin[j];k_idx < G->begin[j+1] ; k_idx ++) {
            node_t s = G_member[j];
	    node_t k = G->node_idx [k_idx];
	    t += (s-(k+2));
	    t = t--;
          }
	}
#pragma omp atomic
      tShared += t;
    }
    tShared /= 1000;
  }

  printf("TSHared %d\n", tShared);
  endtracking();
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  free(G_member);
}






#define numTimes 7

int runalgo(int argc,char** argv) {

  int i;
  //setaffinity();
  graph* G = readGraph(argv[1], argv[2]);
  for(i = 0;i< numTimes; i++) {
    printf("Run %d \n", i);
    train7(G,i);
    sleep(2);
  }
  return 0;
}

inline void kernel(graph *G) {
  
}


