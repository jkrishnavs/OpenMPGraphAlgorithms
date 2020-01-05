#include "graphEnum.h"
#include "print.h"
#include<stddef.h>
#include"cvector.h"




vector* initVector(vector *v, graphDataType t) {
  v = (vector*) malloc (sizeof(vector));
  v->arr = NULL;
  v->capacity = 0;
  v->size = 0;
  v->type = t;
  return v;
}
  
void clearVector(vector *v) {
  v->size = 0;
}

value_t vectorCapacity(vector *v) {
  return v->capacity;
}

value_t vectorSize(vector *v) {
  return v->size;
}


void vectorReserve(vector *v, size_t t) {
  if((value_t)t < v->capacity)
    return;
  if(v->arr == NULL) {
    v->arr = (value_t*) malloc(sizeof(value_t) * t);
  } else {
    value_t* op = v->arr;
    v->arr = (value_t*) realloc(op, sizeof(value_t) * t);
  }
  assert(v->arr != NULL);
  v->capacity = (value_t) t;
}


inline void push_back(vector *v, value_t val) {
  v->arr[v->size] = val;
  v->size ++;
}

value_t getVectorData(vector *v, value_t pos) {
  if(pos > v->size) {
    printError(ARRAY_ACCESS_OUT_OF_BOUNDS,0,NULL);
    return NIL_VAL;
  }
  return v->arr[pos];
}

void setVectorData(vector *v, value_t pos, value_t val) {
  if(pos > v->capacity) {
    printError(ARRAY_ACCESS_OUT_OF_BOUNDS,0,NULL);
    return;
  }
  v->arr[pos] = val;
  if(pos >= v->size)
    v->size = pos +1;
}
