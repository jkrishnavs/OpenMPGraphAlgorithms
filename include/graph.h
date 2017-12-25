#ifndef GRAPH_HEADER_H
#define GRAPH_HEADER_H

//typedefs


#include<stdbool.h>

#ifdef GM_NODE64
typedef int64_t edge_t;
typedef int64_t node_t;
typedef int64_t value_t;
#define GM_SIZE_CHECK_VAR link_error_becuase_gm_graph_lib_is_configured_as_node64_edge64_but_the_application_is_not
#else
typedef int32_t value_t;
typedef int32_t edge_t;
typedef int32_t node_t;
#define GM_SIZE_CHECK_VAR link_error_becuase_gm_graph_lib_is_configured_as_node32_edge32_but_the_application_is_not
#endif


typedef edge_t edge_id;
typedef node_t node_id;

struct edge_dest_t  // for flexible graph representation
{
    node_id dest;
    edge_id edge;
};

static const node_t NIL_NODE = (node_t) -1;
static const edge_t NIL_EDGE = (edge_t) -1;
static const value_t NIL_VAL = (value_t) -1;
/*
  Graph is stored in CSR format.
  Graph has N nodes and M edges.
 */

struct graph
{

  edge_t* begin; /* row ptr of the edges. Array size = N. */
  node_t* node_idx; /* sparse column. Each value represents the destination of the edge. Array size = M. */
  node_t* node_idx_src;  /* sparse column. Created for fast processing. Each value represents the source of the edge. Array size = M. */


  /* Below graph  representation of reverse graph (Where all edges in the graph is reversed).
     Some algorithms may want to process the information from desitnation to source.
     example: PageRank.
     Helper function: createReverseEdges(graph * G);
     In such cases we can  make use of this.*/
  edge_t* r_begin; 
  node_t* r_node_idx;
  node_t* r_node_idx_src;


  /** Mapping created after semi sorting. 
      Helper Function: semisorting
   **/
  edge_t* e_idx2idx;
  edge_t* e_revidx2idx;


  node_t numNodes;
  edge_t numEdges;

  bool reverseEdge;
  bool frozen;
  bool directed;
  bool semiSorted;  
};

typdef struct graph graph;

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
  G->e_rev2idx = NULL;
  reverseEdge = false;
  frozen = false;
  directed = false;
  semiSorted = false;
  numNodes = 0;
  numEdges = 0;
  return G;
}


bool hasEdgeTo(node_t src, node_t dest) {
  // TODO
}
edge_t get_edge_idx_for_src_dest(node_t src, node_t dest) {
  // TODO
}
bool has_separate_edge_idx(graph *G) {
  return (G->e_id2idx != NULL);
}

bool is_edge_source_ready(graph * G) {
  return (G->node_idx_src != NULL);
}

bool isNode( graph *G, node_id n) {
  return (n< G->numNodes);
}

bool isEdge(graph *G, edge_id e) {
  return (e < F->numEdges);
}



#endif
