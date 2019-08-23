#ifndef __cplusplus
#define _GNU_SOURCE
#endif
#ifndef PARSE_GRAPH_H
#define PARSE_GRAPH_H
#include<stdio.h>
#include "graph.h"

int isNumber(char str[]);
value_t numLinesfromCur(FILE *f);
void skipNlines(FILE* f, int n);
void writeBackGraph(graph *G, const char* filename);
void writeEdgeListGraph(graph *G, const char* filename);
graph* parseGraph(const char* filename, const char* graphformat);
graph* parseedgeListGraph(const char* filename, const char* graphformat);


#endif
