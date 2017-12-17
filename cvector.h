#ifndef C_VECTOR_H
#define C_VECTOR_H
/**
   A vector can be required for either 
 **/

#include "graphEnum.h"
#include "print.h"

typedef struct vector{
  value_t* arr;
  value_t capacity;
  value_t size;
  graphDataType type; 
} vector;


void initVector(vector &v, graphDataType t) {
  v.arr = NULL;
  v.capacity = 0;
  v.size = 0;
  v.type = type;
}
  
inline void clearVector(vector &v) {
  v.size = 0;
}

inline value_t vectorCapacity(vector &v) {
  return v.capaity;
}

inline value_t vectorSize(vector &v) {
  return v.size;
}

void vectorReserve(vector &v, size_t t) {
  if(v.capacity < (value_t)t)
    return;
  if(v.arr == NULL) {
    v.arr = (value_t*) malloc(sizeof(value_t) * t);
  } else {
    v.arr = (value_t*) realloc(sizeof(value_t) * t);
  }
  assert(v.arr != NULL);
  v.capacity = (value_t) t;
}

inline void push_back(vector &v, value_t val) {
  v.arr[v.size] = val;
  v.size ++;
}

inline value_t data(vector &v, value_t pos) {
  if(pos > v.size) {
    printError(ARRAY_ACCESS_OUT_OF_BOUNDS);
    return NIL_VAL;
  }
  return v.arr[pos];
}
#endif
