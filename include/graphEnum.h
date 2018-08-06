#ifndef _GRAPH_ENUM_H
#define _GRAPH_ENUM_H

typedef enum executionSection {
  GRAPHREAD,
  REVERSE_EDGE_CREATION,
  GRAPHWRITE,
  ALGO_KERNEL,
  OVERALL
}executionSection;


typedef enum errorCodes {
  GRAPH_FILE_NOT_FOUND,
  CONFIG_FILE_NOT_FOUND,
  GRAPH_FILE_NOT_CREATED,
  OUT_OF_MEMORY,
  ARRAY_ACCESS_OUT_OF_BOUNDS,
  TASKLOOP_NOTENABLED,
  INCORRECT_ARG_LIST
}errorCodes;


typedef enum  graphDataType{
  NODE_T,
  EDGE_T,
  INT32,
  INT64,
  
}graphDataType;



#endif
