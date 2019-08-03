#ifndef GRAPH_HEADER_H
#define GRAPH_HEADER_H

//typedefs


#include<stdbool.h>
#include<stdint.h>
#include<stdlib.h>
#include<stddef.h>
#include<assert.h>

#ifdef GM_NODE64
typedef int edge_t;
typedef int node_t;
typedef int value_t;
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

  edge_t* begin; /* row ptr of the edges. Array size = N +1. */
  node_t* node_idx; /* sparse column. Each value represents the destination of the edge. Array size = M. */
  node_t* node_idx_src;  /* sparse column. Created for fast processing. Each value represents the source of the edge. Array size = M. */

  int* weights; /* Edge weights */
  

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
  bool weighted;
};

typedef struct graph graph;

void deleteGraph(graph *G);
graph* createGraph();
bool hasEdgeTo(graph *G, node_t src, node_t dest);
edge_t get_edge_idx_for_src_dest(graph *G, node_t src, node_t dest);
bool has_separate_edge_idx(graph *G);
bool is_edge_source_ready(graph * G);

bool isNode( graph *G, node_id n);
bool isEdge(graph *G, edge_id e);


#endif
