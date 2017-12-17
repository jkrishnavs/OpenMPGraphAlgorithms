#ifndef _GRAPH_ENUM_H
#define _GRAPH_ENUM_H

enum executionSection {
  GRAPHREAD,
  REVERSE_EGDE_CREATION,
  GRAPHWRITE,
  ALGOKERNEL,
  OVERALL
};


enum errorCodes {
  GRAPH_FILE_NOT_FOUND,
  OUT_OF_MEMORY,
};



#endif
