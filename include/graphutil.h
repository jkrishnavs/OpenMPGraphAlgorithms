#ifndef GRAPH_UTIL_H
#define GRAPH_UTIL_H

#include "graph.h"
#include "atomics.h"
#include "print.h"
#include "vector.h"

/*****
 * These are set of utility functions which follows the same 
 * algorithm as in GreenMarl code. Some of these functions were
 * part of the base gm_graph class. For more details 
 * see original functions in GreenMarl code
 *****/

/**
 * see gm_graph::create_reverse_edges()...
 */ 
void createReverseEdges(graph* G) {

  
  if(G->reverseEdge == true) return;
  if(G->frozen == false) freezeGraph(G);
  node_t nodes = G->numNodes;
  node_t edges = G->numEdges;
  G->r_begin = (edge_t*) malloc ((nodes+1)* sizeof(edge_t));
  assert(G->r_begin != NULL);
    
  G->r_node_idx = (node_t*) malloc ((edges)* sizeof(node_t));
  assert(G->r_node_idx != NULL);
  edge_t* loc = (egde_t*) malloc ((nodes)* sizeof(edge_t));
  assert(loc != NULL);
  // initialize
  
#pragma omp parallel for
  for (node_t i = 0; i < nodes; i++) {
    G->r_begin[i] = 0;
  }


#pragma omp parallel for schedule(dynamic,128)
  for (node_t i = 0; i < nodes; i++) {
    for (edge_t e = G->begin[i]; e < G->begin[i + 1]; e++) {
      node_t dest = G->node_idx[e];
      edge_t location = ATOMIC_INC(&(G->r_begin[dest]));
      loc[e] = location;	
    }
  }
#if !USE_PARALLEL_PREFIXSUM
  edge_t sum = 0;
  for (node_t i = 0; i < nodes; i++) {
    edge_t sum_old = sum;
    sum = sum + G->r_begin[i];
    G->r_begin[i] = sum_old;
  }
  G->r_begin[nodes] = sum;
#else
  // TODO
#endif
  
  
    // now update destination
#pragma omp parallel for schedule(dynamic,128)
  for (node_t i = 0; i < nodes; i++) {
    for (edge_t e = G->begin[i]; e < G->begin[i + 1]; e++) {
      node_t dest = G->node_idx[e];
      edge_t r_edge_idx = G->r_begin[dest] + loc[e];
      G->r_node_idx[r_edge_idx] = i;
    }
  }

  
  
  G->reverseEdge = true;  
  if (G->semiSorted) semisortReverse(G);
  if (G->node_idx_src != NULL) prepareEdgeSourcReverse();

  free(loc);     
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

  semisortMain(G->numNodes, G->numEdges, G->begin, G->node_idx, G->e_idx2idx, G->e_idx2id);
  
#pragma omp parallel for
  for (edge_t j = 0; j < G->numEdges; j++) {
    edge_t id = G->e_idx2id[j];
    G->e_id2idx[id] = j;
  }
  

  if (G->reverseEdge) {
    do_semi_sort_reverse(G);
  }

  G->semiSorted = true;
}

/**
 * see gm_graph::do_semi_sort_reverse()...
 */
void semisortReverse(graph *G) {
  assert(G->semiSorted == true);
  semisortMain(G->numNodes, G->numEdges, G->r_begin, G->r_node_idx, G->e_rev2idx, NULL);
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
    
    vector index;
    initVector(index, NODE_T);
    vector destCopy;
    initVector(destCopy, NODE_T);
    vector auxCopy;
    initVector(auxCopy, EDGE_T);
    vector aux2Copy;
    initVector(aux2Copy, EDGE_T); 

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
	vectorReserve(index);
	vectorReserve(destCopy);
	vectorReserve(auxCopy);
	vectorReserve(aux2Copy);
	
      }
      for(edge_t j=0;j < sz; j++) {
	vectorData(index, j,j);
	vectorData(destCopy, j, destLocal[j]);
	vectorData(auxCopy, j, auxLocal[j]);
	if (aux2 != NULL)
	  vectorData(aux2Copy, j, aux2Local[j]);
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
	destLocal[j] = vectorData(destCopy, vectorData(index, j));
	auxLocal[j] = vectorData(auxCopy, vectorData(index, j) );
	if (aux2 != NULL)
	  aux2Local[j] = vectorData(aux2Copy, vectorData(index,j));
      }
    }
  }
}

#endif
