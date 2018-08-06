#ifndef MAIN_FUNCTION_H
#define MAIN_FUNCTION_H
#include "graph.h"
#include "parsegraph.h"
#include "graphutil.h"
#include<sys/time.h>


int runalgo(int argc,char** argv);
void kernel(graph *G);

graph* readGraph(const char* filename, const char* graphformat) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  graph* G = parseedgeListGraph(filename, graphformat); 
  gettimeofday(&end, NULL);
  printTiming(GRAPHREAD,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  return G;
}


void createReverseEdge(graph* G) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  createReverseEdges(G);
  gettimeofday(&end, NULL);
  printTiming(REVERSE_EDGE_CREATION,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void writeGraph(graph* G, const char* filename) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  writeBackGraph(G, filename);
  gettimeofday(&end, NULL);
  printTiming(GRAPHWRITE,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));

}


/***
  * The main Function.
 **/
int main(int argc, char** args) {
  struct timeval start, end;
  gettimeofday(&start, NULL);
  //omp_set_num_threads(4);	/*  */
#pragma omp parallel
  {

#pragma omp master
    {
      printf("The number of threads is %d \n", omp_get_num_threads());
    }
  }
  runalgo(argc, args);		/*  */
  gettimeofday(&end, NULL);
  printTiming(OVERALL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


void runKernel(graph* G) {
 struct timeval start, end;
  gettimeofday(&start, NULL);
  kernel(G);
  gettimeofday(&end, NULL);
  printTiming(ALGO_KERNEL,((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
}


#endif
