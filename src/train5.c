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
  
}


void trainaed(graph* G, int id) {
  printf("The train aed %d \n", id);
  char title[50];
  sprintf(title, "trainaed_%d.csv",id);
  gettimeofday(&start, NULL);
  inittracking(title);
  double* aed = (double*) malloc (G->numNodes * sizeof(double));

#pragma omp parallel for
  for (node_t t = 0; t < G->numNodes; t ++) 
  {
    aed[t] = t/ 1024;
  } 
  for(int abc=0; abc < REPEAT; abc ++) { 
#pragma omp parallel
    {
      double d = 0.0;
      int final = 0;
#pragma omp for schedule(dynamic,1024)
      for (node_t t = 0; t < G->numNodes; t ++) 
	{
	  for (edge_t w_idx = G->begin[t];w_idx < G->begin[t+1] ; w_idx ++) 
	    {
	      node_t w = G->node_idx [w_idx];
	      d += aed[w];
	      final += (G->begin[w+1] - G->begin[w]);
	    }
	  final = final % G->numNodes;
	}
#pragma omp atomic write
      aed[final] = d;
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
    trainaed(G,i);
    sleep(2);
  }
  return 0;
}

inline void kernel(graph *G) {
  
}


