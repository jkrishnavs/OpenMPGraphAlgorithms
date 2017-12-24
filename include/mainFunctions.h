#ifndef MAIN_FUNCTION_H
#define MAIN_FUNCTION_H
#include "graph.h"
#include "parseGraph.h"
#include "graphutil.h"
#include<time.h>
graph* readGraph(static const char[] filename) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  graph* G = parseGraph(static const char[] filename); 
  gettimeofday(&end, NULL);
  printtiming(GRAPHREAD,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void createReverseEdge(graph* G) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  createReverseEdges(G);
  gettimeofday(&end, NULL);
  printtiming(REVERSE_EDGE_CREATION,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void writeGraph(graph* G, static const char[] filename) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  writeBackGraph(G, filename);
  gettimeofday(&end, NULL);
  printtiming(GRAPHWRITE,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));

}


/***
  * The main Function.
 **/
int main(int argc, char** args) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  runalgo(argc, args);
  gettimeofday(&end, NULL);
  printtiming(OVERALL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void runKernel(int argc, char** args) {
 struct timeval start, end;
  gettimeofday(&start, NULL);
  *fn(argc, args);
  gettimeofday(&end, NULL);
  printtiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


#endif
