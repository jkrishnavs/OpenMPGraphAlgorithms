#ifndef GRAPH_UTIL_H
#define GRAPH_UTIL_H

#include "graph.h"
#include "atomics.h"
#include "print.h"
#include "cvector.h"


/*****
 * These are set of utility functions which follows the same 
 * algorithm as in GreenMarl code. Some of these functions were
 * part of the base gm_graph class. For more details 
 * see original functions in GreenMarl code
 *****/

/* Declarations */
bool isNeighbour(graph *G, node_t src, node_t dest);
void semisort(graph *G);
void semisortReverse(graph *G);
void createReverseEdges(graph* G);
void prepareEdgeSource(graph *G);
void prepareEdgeSourceReverse(graph *G);
void freezeGraph(graph * G);
void semisortMain(node_t N, edge_t M, edge_t* begin, node_t* dest, edge_t* aux, edge_t* aux2);

bool isNeighbour(graph *G, node_t src, node_t dest) {
  bool neighbour = false;
  for(edge_t s = G->begin[src]; s < G->begin[src+1]; s++) {
    node_t y = G->node_idx[s];
    if(y == dest) {
      neighbour = true;
      break;
    }
  }
  return neighbour;
}

/**
 * see gm_graph::do_semi_sort()...
 **/
void semisort(graph *G) {
  if(G->semiSorted == true) return;
  if(G->frozen == false) freezeGraph(G);


  G->e_idx2idx = (edge_t*) malloc(sizeof(edge_t)* G->numEdges);
  assert(G->e_idx2idx != NULL);
#pragma omp parallel for schedule(dynamic,128)
  for (node_t i = 0; i < G->numNodes; i++) {
    for (edge_t j = G->begin[i]; j < G->begin[i + 1]; j++) {
      G->e_idx2idx[j] = j;    
    }
  }

  semisortMain(G->numNodes, G->numEdges, G->begin, G->node_idx, G->e_idx2idx, G->e_idx2idx);

  /* TODO uncomment this region 
   if we require e_ixd2idx later
  */
/* #pragma omp parallel for */
/*   for (edge_t j = 0; j < G->numEdges; j++) { */
/*     edge_t id = G->e_idx2idx[j]; */
/*     G->e_idx2idx[id] = j; */
/*   } */

  /* TODO: .. and comment this */
  free(G->e_idx2idx);
  G->e_idx2idx = NULL;

  
  if (G->reverseEdge) {
    semisortReverse(G);
  }

  G->semiSorted = true;
}

/**
 * see gm_graph::do_semi_sort_reverse()...
 */
void semisortReverse(graph *G) {
  assert(G->semiSorted == true);

  if(G->e_revidx2idx == NULL) {
    G->e_revidx2idx = (edge_t*) malloc(sizeof(edge_t)* G->numEdges);
    assert(G->e_revidx2idx != NULL);
#pragma omp parallel for schedule(dynamic,128)
    for (node_t i = 0; i < G->numNodes; i++) {
      for (edge_t j = G->r_begin[i]; j < G->r_begin[i + 1]; j++) {
	G->e_revidx2idx[j] = j;    
      }
    }
  }
  
  semisortMain(G->numNodes, G->numEdges, G->r_begin, G->r_node_idx, G->e_revidx2idx, NULL);
  /* TODO : If we are using  e_revidx2idx uncomment the following region*/

/* #pragma omp parallel for schedule(dynamic,128) */
/*   for (node_t i = 0; i < G->numNodes; i++) { */
/*     for (edge_t j = G->r_begin[i]; j < G->r_begin[i + 1]; j++) { */
/*       G->e_revidx2idx[j] = j;     */
/*     } */
/*   } */

  /* TODO: .. and comment this region */
  free(G->e_revidx2idx);
  G->e_revidx2idx = NULL;
  
}




/**
 * see gm_graph::create_reverse_edges()...
 */ 
void createReverseEdges(graph* G) {

  
  if(G->reverseEdge == true) return;
  if(G->frozen == false) freezeGraph(G);
  node_t nodes = G->numNodes;
  edge_t edges = G->numEdges;
  G->r_begin = (edge_t*) malloc ((nodes+1)* sizeof(edge_t));
  assert(G->r_begin != NULL);
    
  G->r_node_idx = (node_t*) malloc ((edges)* sizeof(node_t));
  assert(G->r_node_idx != NULL);
  
  /* edge_t* loc = (edge_t*) malloc ((nodes)* sizeof(edge_t)); */
  /* assert(loc != NULL); */
  // initialize

  node_t* relLoc = (node_t*) malloc ((nodes)* sizeof(node_t));

  node_t i;
#pragma omp parallel for
  for (i = 0; i < nodes; i++) {
    G->r_begin[i] = 0;
    relLoc[i] = 0;
  }
  
#pragma omp parallel for schedule(dynamic,128)
  for (i = 0; i < nodes; i++) {
    for (edge_t e = G->begin[i]; e < G->begin[i + 1]; e++) {
      node_t dest = G->node_idx[e];
      #pragma omp atomic
      G->r_begin[dest]++;
    }
  }
  
#if !USE_PARALLEL_PREFIXSUM
  for (i = 0; i <= nodes; i++) {
    G->r_begin[i] += G->r_begin[i-1];
  }
#else
  // TODO
#endif


  
  /* Our Experiments show that atleast for a smaller 
     number of hardware threads a sequential algorithm 
     is efficient in creating the reverse edges (no synchronization costs). 
     The reverse edges created using sequential set will
     be inadvertantly sorted. */


  
  for (i = 0; i < nodes; i++) {
    edge_t e;
    for (e = G->begin[i]; e < G->begin[i + 1]; e++) {
      node_t dest = G->node_idx[e];
      edge_t r_edge_idx = G->r_begin[dest]+ (relLoc[dest]++);
      G->r_node_idx[r_edge_idx] = i;
    }
  }

  
  
  G->reverseEdge = true;  
  // if (G->semiSorted) semisortReverse(G);
  // if (G->node_idx_src != NULL) prepareEdgeSourceReverse(G);

  free(relLoc);     
}




/**
 * see gm_graph::prepare_edge_source()...
 */
void prepareEdgeSource(graph *G) {
  if(G->node_idx_src != NULL) return;
  
  G->node_idx_src = (node_t*) malloc(G->numEdges*sizeof(node_t));
  assert(G->node_idx_src != NULL);

#pragma omp parallel for schedule(dynamic,128)
  for (node_t i = 0; i < G->numNodes; i++) {
    for (edge_t j = G->begin[i]; j < G->begin[i + 1]; j++) {
      G->node_idx_src[j] = i;
    }
  }
  if(G->reverseEdge ==  true)
    prepareEdgeSourceReverse(G);
}

/**
 * see gm_graph::prepare_edge_source_reverse()...
 */
void prepareEdgeSourceReverse(graph *G) {
  assert(G->node_idx_src != NULL);
  G->r_node_idx_src = (node_t*) malloc(G->numEdges*sizeof(node_t));
  assert(G->r_node_idx_src != NULL);
  
#pragma omp parallel for schedule(dynamic,128)
  for (node_t i = 0; i < G->numNodes; i++) {
    for (edge_t j = G->r_begin[i]; j < G->r_begin[i + 1]; j++) {
      G->r_node_idx_src[j] = i;
    }
  }
}

/**
 * see gm_graph::freeze()...
 */
void freezeGraph(graph * G) {
  if(G->frozen == true) return;  
  /***
   * TODO :We are not handling flexible graphs for time being.
   **/
  
  G->frozen = true;
  G->semiSorted = false;
  G->reverseEdge = false;
  semisort(G);
#ifdef NON_NORMALIZED_GRAPH
  // make sure the node_idx_src is filled properly.
#endif
}


/**
 * see semi_sort_main()...
 */
void semisortMain(node_t N, edge_t M, edge_t* begin, node_t* dest, edge_t* aux, edge_t* aux2) {
#pragma omp parallel
  {
    
    vector* index = NULL;
    index = initVector(index, NODE_T);
    vector* destCopy = NULL;
    destCopy = initVector(destCopy, NODE_T);
    vector* auxCopy = NULL;
    auxCopy = initVector(auxCopy, EDGE_T);
    vector* aux2Copy = NULL;
    aux2Copy = initVector(aux2Copy, EDGE_T); 

#pragma omp for schedule(dynamic,4096) nowait
    for (node_t i = 0; i < N; i++) {
      clearVector(index);
      clearVector(destCopy);
      clearVector(auxCopy);
      clearVector(aux2Copy);
      edge_t sz = begin[i+1] - begin[i];
      node_t* destLocal = dest + begin[i];
      edge_t* auxLocal =  aux + begin[i];
      edge_t* aux2Local = (aux2 == NULL)?NULL:aux2 + begin[i];
      
      if (vectorCapacity(index) < (size_t) sz) {
	vectorReserve(index, sz);
	vectorReserve(destCopy, sz);
	vectorReserve(auxCopy, sz);
	vectorReserve(aux2Copy, sz);
	
      }
      for(edge_t j=0;j < sz; j++) {
	setVectorData(index, j,j);
	setVectorData(destCopy, j, destLocal[j]);
	setVectorData(auxCopy, j, auxLocal[j]); /*  */
	if (aux2 != NULL)
	  setVectorData(aux2Copy, j, aux2Local[j]);
      }
      // TODO sort.
      // sort indicies
      /***
       * The stl sort of C++ (used by GreenMarl) is faster than 
       * the C quick sort.
       * much faster. Should we go for hand tuned sort here? 
       * Details: http://www.geeksforgeeks.org/c-qsort-vs-c-sort/
       * std::sort(index.data(), index.data()+sz, _gm_sort_indices<node_t>(dest_copy.data()));
       **/

      for(edge_t j=0;j < sz; j++) {
	destLocal[j] = getVectorData(destCopy, getVectorData(index, j));
	auxLocal[j] = getVectorData(auxCopy, getVectorData(index, j) );
	if (aux2 != NULL)
	  aux2Local[j] = getVectorData(aux2Copy, getVectorData(index,j));
      }
    }
  }
}

#endif
