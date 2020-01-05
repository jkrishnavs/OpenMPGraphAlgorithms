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
// This is a very large kernel, so repeat only once.
#define REPEAT 1

long long iters[8];

struct timeval start, end;

// We define all additional paramemter here
void setaffinity() {
  
}


void train11(graph* G, int id) {
  printf("The train 11 %d \n", id);
  char title[50];
  sprintf(title, "train11_%d.csv",id);
  gettimeofday(&start, NULL);
  inittracking(title);
  double* aed = (double*) malloc (G->numNodes * sizeof(double));
  int T = -10000;
#pragma omp parallel for
  for (node_t t = 0; t < G->numNodes; t ++) 
  {
    aed[t] = t/ 1024;
  } 
  for(int abc=0; abc < REPEAT; abc ++) { 
#pragma omp parallel
    {
      int T_prv = 0;
      node_t v;
#pragma omp for
      for (v = 0; v < G->numNodes; v ++) {
	edge_t u_idx;
	for (u_idx = G->begin[v]; u_idx < G->begin[v+1]; u_idx ++) {
	  node_t u = G->node_idx [u_idx];
	  if (u > v) {
	    edge_t w_idx;
	    for (w_idx = G->begin[u]; w_idx < G->begin[u+1]; u_idx ++) {
	      node_t w = G->node_idx [w_idx];
	      if (w > u) {
		edge_t s;
		for(s= G->begin[w]; s < G->begin[w+1]; s++) {
		  node_t y = G->node_idx[s];
		  if(y == u) {
		    T_prv = T_prv + 1 ;
		    break;
		  }
		}
	      }
	    }
	  }
	}
      }
#pragma omp atomic
      T += T_prv;
    }
  }
  for(int i=0; i< 10;i++) {
    printf("aed[%d] =  %f \n", i,aed[i]);
  }
  free(aed);
  endtracking();
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));

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
    train11(G,i);
    sleep(2);
  }
  return 0;
}

inline void kernel(graph *G) {
  
}


