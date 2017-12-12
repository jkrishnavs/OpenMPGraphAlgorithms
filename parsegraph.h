#include "graph.h"
#include<string.h>
#include<file.h>
#include <ctype.h>


int isNumber(char[] str) {
  i=0;
  while(s[i] != '\0') {
    if (isdigit(s[i]) == 0) return 0;
  }
  return 1;
}

value_t numLinesfromCur(FILE *f) {
  int lines  = 0;
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
  while(!feof(f) && counter < n) {
    fscanf(f, "%*[^\n]\n", NULL);
    counter++;
  }
}


graph* parseGraph(static const char[] filename) {
  FILE *f;

  f = fopen(filename, "r");
  if(f == NULL) {
    printf("Error: File Not Found \n");
    return NULL;
  }

  
  /**
     TODO
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
  int isNumber = 0;
  // Parse all nodes
  
  while(f != EOF && isnumber == 0) {
    fscanf(f,"%d %s", x, str);
    isNumber = isNumber(str);
    if(isNumber == 0) {
      /***
	  TODO: We need to revist this.
	  Now we assume all inpt nodes contains 0 to N-1 nodes.
       **/
      assert(x == counter);
      counter++;
    }  
  }
  G->numNodes = counter;

  G->begin = (edge_t*) malloc (sizeof(edge_t) * counter);
  assert(G->begin != NULL);
  /**get the number of edges **/
  G->numEdges = (edge_t) numLinesfromCur(f) + 1;


  fseek(f,0,SEEK_SET); // rewind
  skipNlines(f, counter);  
  assert(f != EOF);

  G->node_idx = (node_t*) malloc(sizeof(node_t) * G->num_Edges);

  assert(G->node_idx != NULL);


  edge_t edgeid = 0;
  node_t y;
  node_t curSource = 0;
  G->begin[curSource] = 0;
  while(f != EOF) {   
    fscanf(f, "%d %d", x,y);
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

  fclose(f);

  return G;
}


