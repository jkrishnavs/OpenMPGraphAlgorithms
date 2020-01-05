#ifndef GRAPH_PROP_H
#define GRAPH_PROP_H
typedef struct graph graph;

double avgEdgeDistance(graph *G);
double avgClusterCoeff(graph *G);
node_t diameter(graph *G);
double sparsity(graph *G);
double sccIndex(graph *G);
double triangle_counting(graph *G);

#endif

