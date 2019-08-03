#ifndef GRAPH_GENERATORS_H
#define GRAPH_GENERATORS_H
#include <sprng_cpp.h>
#include "graph.h"
#include "graphProperty.hh"


/** C functions **/
graph* allocateMemoryforGraph(node_t n, edge_t m, bool weighted);
void doubleMergeSort(node_t* l1, node_t* l2, edge_t left, edge_t right);
void merge(node_t* l1, node_t* l2, edge_t left, edge_t mid, edge_t right);
void choosePartition(RMATdata& r, node_t* u, node_t* v, node_t step, Sprng* stream1);
void varyParams(double* a, double* b, double* c, double* d, Sprng * stream3, Sprng* stream4);


graph* randomGenerator(GraphProperty& p);
graph* erdosRenyiGenerator(GraphProperty& p);
graph* rmatGenerator(GraphProperty& p);
graph* SSCAGenerator(GraphProperty& p);
//graph* propertyControlledGraphGenerator(GraphProperty& p);

graph* StocasticBlockModel(GraphProperty& p);
graph* lowdegreeNetworkgraph(GraphProperty& p);
graph* diagonalGraphGenerator(bool randomizedmidpoint, GraphProperty& p);
graph* SBMGraphGenerator(double* probabilityVector, double *edgeProbabilityMatrix, node_t k, GraphProperty& p);
graph* symmetricStocasticBlockModel(node_t k, double A, double B, GraphProperty& p);

graph* diagonalGraphGenerator(bool randomizedmidpoint, GraphProperty& p);
graph* adjustAED(graph *G, gdata& data, Sprng* stream, int itrs, int flag);
graph* changeClustering(graph* G, gdata& data, Sprng* stream, int itrs, int flag);
graph* changeDegreesd(graph* G, gdata& data, Sprng* stream, int itrs, int flag);

#endif
