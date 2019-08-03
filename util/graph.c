#include <stdio.h>
#include "graph.h"


void deleteGraph(graph *G) {

  free(G->begin);
  free(G->node_idx);
  if(G->node_idx_src != NULL)
    free(G->node_idx_src);
  if(G->weights != NULL)
    free(G->weights);
  if(G->r_begin != NULL)
    free(G->r_begin);
  if(G->r_node_idx != NULL)
    free(G->r_node_idx);
  if(G->r_node_idx_src != NULL)
    free(G->r_node_idx_src);
  if(G->e_idx2idx != NULL)
    free(G->e_idx2idx);
  if(G->e_revidx2idx != NULL)
    free(G->e_revidx2idx);
  free(G);
}

graph* createGraph() {
  graph* G  = (graph*) malloc(sizeof(graph));
  assert(G != NULL);
  G->begin = NULL;
  G->node_idx = NULL;
  G->node_idx_src = NULL;  
  G->r_begin = NULL;
  G->r_node_idx = NULL;
  G->r_node_idx_src = NULL;
  G->e_idx2idx = NULL;
  G->e_revidx2idx = NULL;
  G->weights = NULL;
  G->reverseEdge = false;
  G->frozen = false;
  G->directed = false;
  G->semiSorted = false;
  G->weighted = false;
  G->numNodes = 0;
  G->numEdges = 0;
  return G;
}


bool hasEdgeTo(graph *G, node_t src, node_t dest) {
  bool hasEdge = false;
  for(edge_t s = G->begin[src]; s < G->begin[src+1]; s++) {
    node_t y = G->node_idx[s];
    if(y == dest) {
      hasEdge = true;
      break;
    }
  }
  return hasEdge;
}

edge_t get_edge_idx_for_src_dest(graph *G, node_t src, node_t dest) {
  edge_t edge = NIL_EDGE;
  for(edge_t s = G->begin[src]; s < G->begin[src+1]; s++) {
    node_t y = G->node_idx[s];
    if(y == dest) {
      edge = s;
      break;
    }
  }
  return edge;
}
bool has_separate_edge_idx(graph *G) {
  return (G->e_idx2idx != NULL);
}

bool is_edge_source_ready(graph * G) {
  return (G->node_idx_src != NULL);
}

bool isNode( graph *G, node_id n) {
  return (n< G->numNodes);
}

bool isEdge(graph *G, edge_id e) {
  return (e < G->numEdges);
}
