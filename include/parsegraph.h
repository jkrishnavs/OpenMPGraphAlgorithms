#ifndef PARSE_GRAPH_H
#define PARSE_GRAPH_H
#include "graph.h"
#include "graphEnum.h"
#include "print.h"
#include <string.h>
#include <stdio.h>
#include <ctype.h>


int isNumber(char str[]) {
  int i=0;
  while(str[i] != '\0') {
    if (isdigit(str[i]) == 0) return 0;
    i++;
  }
  return 1;
}

value_t numLinesfromCur(FILE *f) {
  int lines  = 0;
  char ch;
  while(!feof(f)) {
      ch = fgetc(f);
      if(ch == '\n')
	{
	  lines++;
	}
  }
  return lines;
}

void skipNlines(FILE* f, int n) {
  int counter  = 0;
  int r = 0;
  while(!feof(f) && counter < n) {
    r = fscanf(f,"%*[^\n]\n");
    assert(r!= EOF);
    counter++;
  }
}


void writeBackGraph(graph *G, const char* filename) {
  /* TODO if requied */
}


graph* parseGraph(const char* filename) {
  FILE *f;

  f = fopen(filename, "r");
  if(f == NULL) {
    printError(GRAPH_FILE_NOT_FOUND,0, NULL);
    return NULL;
  }

  
  /**
     TODO if required.
   We assume the edge file format to be listing all 
   nodes initially in format (one in each line).
   nodeid *
   followed by edges in format (one in each line).
   src dest
   We have added some asserts to check if this format is
   indeed followed, else we will comeback to change this.
   **/

  
  graph *G = createGraph();

  node_t x;
  char str[20];

  node_t counter  = 0;
  int isnumber = 0;
  // Parse all nodes
  int r = 1;
  while(r != EOF && isnumber == 0) {
    r = fscanf(f,"%d %s", &x, str);
    isnumber = isNumber(str);
    if(isnumber == 0) {
      /***
	  TODO: if required 
	  We need to revisit this if assert fails.
	  Now we assume all input nodes contains 0 to N-1 nodes.
       **/
      assert(x == counter);
      counter++;
    }  
  }
  G->numNodes = counter;

  G->begin = (edge_t*) malloc (sizeof(edge_t) * (G->numNodes+1));
  assert(G->begin != NULL);
  /**get the number of edges **/
  G->numEdges = (edge_t) numLinesfromCur(f);
  printf("The number of Nodes is %d\n", G->numNodes);
  printf("The number of Edges is %d\n", G->numEdges);

  fseek(f,0,SEEK_SET); // rewind
  skipNlines(f, counter);  
  assert(!feof(f));

  G->node_idx = (node_t*) malloc(sizeof(node_t) * G->numEdges);

  assert(G->node_idx != NULL);


  r = 1;

  /**
   * TODO if required: 
   * We currently assume that all edges are ordered.
   * assert (curSource < x) will fail.
   **/ 

  edge_t edgeid = 0;
  node_t y;
  node_t curSource = 0;
  G->begin[curSource] = 0;
 
  while(r != EOF) {   
    r = fscanf(f, "%d %d", &x, &y);
    if(r != EOF) {
      //printf("The edge id is %d from %d to %d \n",edgeid,x,y);
      G->node_idx[edgeid] = y;
      // format check
      if(x != curSource) {
	assert(curSource < x);
 	while(x != curSource) {
	  curSource++;
	  G->begin[curSource] = edgeid;
	}
      }
      edgeid++;
    }
  }

  while(curSource < G->numNodes) {
    curSource++;
    G->begin[curSource] = edgeid;
  }

  assert(edgeid == G->numEdges);

  //  printf("The current source is %d\n",curSource);

  fclose(f);

  return G;
}


#endif
