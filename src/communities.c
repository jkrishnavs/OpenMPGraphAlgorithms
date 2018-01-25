/**
 * Original base algorithm from GreenMarl.
 * Otimized and edited by
 * edited by Jyothi Krishna V S.
 */

#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "nodeIntMap.h"
#include "communities.h"

#define NO_OF_ARGS 2


void communities(graph* G);




void output(graph *G) {
  // print output.
  node_t commList[10];
  int commCount[10];
  int i;
  for(i=0;i<10; i++) {
    commList[i] = NIL_NODE;
    commCount[i] = 0;
  }
  int found;
  int curIndex = 0;
  int t;
  for (t = 0; t < G->numNodes; t++) {
    found  = -1;
    for(i = 0; i<10;i++) {
      if(comm[t] == commList[i]) {
	found = i;
	break;
      }
    }
    if(found != -1) {
      commCount[found]++;
    } else if(curIndex < 10) {
      commCount[curIndex] = 1;
      commList[curIndex] = comm[t];
      curIndex++;
    }    
  }
  printf("Community\t#Nodes\t\t(Showing max 10 entries)\n");
  for (i=0; i<10; i++) {
    if(commList[i] != NIL_NODE)
      printf("%d\t\t%d\n", commList[i], commCount[i]);
  }
  free(comm);
  comm = NULL;
}


/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  int flag = 0;
  maxItrs = 10;

  if(argc > 2) {
    maxItrs = atoi(argv[2]);
    if(maxItrs <=0) flag = 1;
  }

  graph* G = readGraph(argv[1]);

  if(G == NULL) {
    flag = 2;
  }
  
  if(flag != 0 ) {
    const char* argList[NO_OF_ARGS] = {" <inputfile> " , "[max_iterations = 10]"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  
  comm = (node_t*) malloc (G->numNodes * sizeof(node_t));
  assert(comm != NULL);
  runKernel(G);
  output(G);
  return 0;
}










inline void kernel(graph *G) {
  communities(G);
}










