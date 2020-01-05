#ifndef C_VECTOR_H
#define C_VECTOR_H
/**
   A vector can be required for either 
 **/

#include"graphEnum.h"
#include"graph.h"

typedef struct vector{
  value_t* arr;
  value_t capacity;
  value_t size;
  graphDataType type; 
} vector;

vector* initVector(vector *v, graphDataType t);
void clearVector(vector *v);
value_t vectorCapacity(vector *v);
value_t vectorSize(vector *v);
void vectorReserve(vector *v, size_t t);
void push_back(vector *v, value_t val);
void setVectorData(vector *v, value_t pos, value_t val);
value_t getVectorData(vector *v, value_t pos);

#endif
