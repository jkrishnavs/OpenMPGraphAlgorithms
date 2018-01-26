/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "powerperformacetracking.h"
#include "print.h"
#include "communities.h"
#include "conduct.h"
#include "pageRank.h"
#include "sssp.h"
#include "triangleCounting.h"
#include <stdlib.h>

#include<unistd.h>

#define NO_OF_ARGS 1

// We define all additional paramemter here

int ssspMaxLength;
int ssspSeed;
void setAdditionalArguments() {
  /* communities */
  maxItrs = 10;
  /* conduct */
  /* pageRank */
  maxIters = 100;
  e = 0.01;
  d = 0.95;
  /* sssp  */
  root =0;
  ssspMaxLength = 100;
  ssspSeed = 0;
  /* triangleCounting */
}



void runcommunities(graph *G) {
  initCommunities();
  struct timeval start, end;
  gettimeofday(&start, NULL);
  communities(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  free(comm);
  comm = NULL;  
}

void runconduct(graph *G) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  conduct(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}

void runpageRank(graph *G) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  runpageRank(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void runtriangleCounting(graph *G) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  runtriangleCounting(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  
}


void runsssp(graph *G) {
  len = (uint32_t*) malloc (G->numEdges * sizeof(uint32_t));
  /*
    Random generation of edge lengths
   */
  srand(ssspSeed);
  for(edge_t e = 0; e < G->numEdges; e++) {
    len[e] = (rand()%ssspMaxLength + 1);
  }
  dist = (uint32_t*) malloc (G->numNodes * sizeof (uint32_t));
  struct timeval start, end;
  gettimeofday(&start, NULL);
  runsssp(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  
  free(len);
  free(dist);
  len = NULL;
  dist = NULL;
}



void output(graph *G) {

  int counter  = 0;
  
  printf("The shortest path to the following nodes from the root Node %d\n", root);
  for (node_t n = 0; n < G->numNodes && counter < 5; n++) {
    if(dist[n] != UINT_MAX && n != root) {
      printf(" %d -> %d : %ud \n", root, n, dist[n]);
      counter ++;
    }
  }  
}



#define numTimes 5

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {

  assert(argc >=2);
  graph* G = readGraph(argv[1]);
  setAdditionalArguments();
  

  dist = (uint32_t*) malloc (G->numNodes * sizeof (uint32_t));
  assert(dist != NULL);
  
  createReverseEdge(G);
  int i;
  for(i = 0;i< numTimes; i++) {
    printf("Run %d \n", i);
    runcommunities(G);
    sleep(30);
    runconduct(G);
    sleep(30);
    runpageRank(G);
    sleep(30);
    runsssp(G);
    sleep(30);
    runtriangleCounting(G);
    sleep(30);
    /* rename the prof files */
    char comm[50];
    sprintf(comm, "communities.%d.csv", i);
    rename("communities.csv", comm);

    char cond[50];
    sprintf(cond, "conduct.%d.csv", i);
    rename("conduct.csv", cond);

    char page[50];
    sprintf(page, "pageRank.%d.csv", i);
    rename("pageRank.csv", page);

    char sssp[50];
    sprintf(comm, "sssp.%d.csv", i);
    rename("sssp.csv",sssp);

    char tri[50];
    sprintf(tri, "triangleCounting.%d.csv", i);
    rename("triangleCounting.csv", tri);
    sleep(30);
    
  }
  return 0;
}



void runAllKernels(graph* G) {
  

}


inline void kernel(graph *G) {
}




