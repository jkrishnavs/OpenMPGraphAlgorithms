#ifndef _GRAPH_ENUM_H
#define _GRAPH_ENUM_H

enum executionSection {
  GRAPHREAD,
  REVERSE_EGDE_CREATION,
  GRAPHWRITE,
  ALGO_KERNEL,
  OVERALL
};


enum errorCodes {
  GRAPH_FILE_NOT_FOUND,
  OUT_OF_MEMORY,
  ARRAY_ACCESS_OUT_OF_BOUNDS,
};


enum  graphDataType{
  NODE_T,
  EDGE_T,
  INT32,
  INT64,
  
};



#endif
