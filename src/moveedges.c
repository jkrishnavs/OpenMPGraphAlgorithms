#include<stdlib.h>
#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
//#include "powerperformacetracking.h"
//#include "communities.h"
#include "graphprop.h"
#include "nodeIntMap.h"
#include <string.h>


#define DEBUG_ON



int adj = 0; // 1 for reverse adjlist or 0 in adjlist.
bool skewed  =  false;

node_t* comm = NULL;

int comp (const void * elem1, const void * elem2) 
{
    int f = *((int*)elem1);
    int s = *((int*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}


graph* moveedges(graph *G) {
  node_t max = G->numNodes;
  unsigned char* flags = (unsigned char*) malloc (sizeof(node_t) * G->numEdges);
  memset(flags, 0, G->numEdges);
  edge_t* nodeT = (edge_t*) malloc (sizeof(edge_t) * G->numNodes);

  for(node_t i =0; i< G->numNodes; i++) {
     nodeT[i]  = 0; //G->begin[i+1] - G->begin[i];
  }
  node_t start, level;
  for(node_t i = 0; i< max;i++) {
    
    do{
      start = rand() % max; // start point to replace
      level = G->begin[start+1] - G->begin[start];
    } while(level ==0);

    node_t end = rand() % level; // end point of edge
    
    if(flags[G->begin[start] + end ] == 0) {
      // new start and end points.
      node_t new_start = rand()% max; 
      //printf("old %d, edge = %d, new = %d \n", start, end, new_start);
      flags[G->begin[start] + end] = 1;
      nodeT[new_start]++;
    } else {
      i--;
    }
  }
  // TODO


  edge_t sub_array_size  =  G->numEdges / G->numNodes;
  node_t* end_points = (node_t*) malloc (sizeof(node_t) * G->numEdges);
  node_t* new_end_points  = (node_t*) malloc (sizeof(node_t) * sub_array_size);

  edge_t ep_index = 0;
  edge_t nep_index = 0;
  edge_t prev_i = 0;
  for(node_t i = 0; i < G->numNodes; i++) {
    if(sub_array_size < nodeT[i]) {
      sub_array_size = nodeT[i];
      new_end_points = (node_t*) realloc (new_end_points,sizeof(node_t) * sub_array_size);
    }
    nep_index = 0;

    // find new end points
    for(edge_t e = 0; e< nodeT[i]; e++) {
      node_t new_end = rand()%max;
      int flag_e = 0;
      for(edge_t e1 = 0;e1 < e; e1++) {
	if(new_end_points[e1] == new_end) {
	  // edge already generated
	  flag_e = 1;
	  break;
	}
      }
      for(node_t old_e = G->begin[i]; old_e < G->begin[i+1]; old_e++) {
	if(G->node_idx[old_e] == new_end && flags[old_e] == 0) {
	  // the edge is already present in original graph
	  flag_e = 1;
	}
      }
      if(new_end == i || flag_e == 1) {
	// self loop , avoid.
	e--;
	continue;
      }
      if(e>=0)
	new_end_points[e] = new_end;
    }
    // sort the end points
    if(nodeT[i] >1)
      qsort(new_end_points, nodeT[i], sizeof(node_t), comp);
    
    // update the array indices
    if(i>0) {
      G->begin[i-1] = prev_i;
      prev_i = ep_index;
    }  

    // place the new nodes in appropriate places.
    for(edge_t ev = G->begin[i]; ev < G->begin[i+1]; ev++) {
      if(!flags[ev]) {
	while(new_end_points[nep_index] < G->node_idx[ev] && nep_index < nodeT[i])
	  end_points[ep_index++] = new_end_points[nep_index++];
	end_points[ep_index++] = G->node_idx[ev];
      }
    } 
    
    // copy the rest of new_end_points;
    while(nep_index < nodeT[i]) {
      end_points[ep_index++] = new_end_points[nep_index++];
    }
    
  }

  G->begin[G->numNodes-1] = prev_i;


  assert(ep_index == G->numEdges);
  node_t* old_edges = G->node_idx;
  G->node_idx  =  end_points;
  free(old_edges);
  free(new_end_points);
  free(flags);
  free(nodeT);
  return G;
}


double avgincident;
double scaledT;
double clusterCoeff;
double aed;
node_t dim;
double sparsityMeasure;
void writeSchema(const char *filename) {
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "Avg Adjecency =  %f \n", avgincident);
  fprintf(fp, "Avg ClusterCoeff = %f \n",clusterCoeff);
  fprintf(fp, "Avg Edge Distance = %f \n", aed);
  fprintf(fp, "Sparsity = %0.7f \n",sparsityMeasure);
  fprintf(fp,"The Scaled percentange Triangles = %f\n\n", scaledT);  
  fclose(fp);
}

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  graph* G = readGraph(argv[1], argv[2]);
  int seed =  0;
  if(argc == 6)
    seed = atoi(argv[5]);
  if(argc < 6) {
    const char* argList[5] = {" <inputfile> " , "graphformat.txt","<outputfile>", "<outputpropfile>", "<seed[default =0]>"};
    printError(INCORRECT_ARG_LIST, 5, argList);
    return -1;

  }
  graph* newG = moveedges(G);

  avgincident = ((double) newG->numEdges )/ newG->numNodes;
  
  writeEdgeListGraph(newG, argv[3]);
  clusterCoeff = avgClusterCoeff(newG);
  aed = avgEdgeDistance(newG);
  dim = diameter(newG);
  sparsityMeasure  = sparsity(newG);
  double T = triangle_counting(newG);

  avgClusterCoeff(newG);
  avgEdgeDistance(newG);
  diameter(newG);
  sparsity(newG);
  triangle_counting(newG);
  scaledT = ((double) T )/ newG->numNodes; 
  sccIndex(newG); 
  writeSchema(argv[4]);
  return 0;
}

inline void kernel(graph *G) {

}
