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

/**** utility functions *****/
graph* readGraph(const char* filename, const char* graphformat);
void createReverseEdge(graph* G);
void writeGraph(graph* G, const char* filename);


/* Declarations */
bool isNeighbour(graph *G, node_t src, node_t dest);
void semisort(graph *G);
void semisortReverse(graph *G);
void createReverseEdges(graph* G);
void prepareEdgeSource(graph *G);
void prepareEdgeSourceReverse(graph *G);
void freezeGraph(graph * G);
void semisortMain(node_t N, edge_t M, edge_t* begin, node_t* dest, edge_t* aux, edge_t* aux2);


#endif
