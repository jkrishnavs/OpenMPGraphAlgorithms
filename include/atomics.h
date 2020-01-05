#ifndef ATOMICS_H
#define ATOMICS_H

#include<omp.h>
#include"graph.h"

inline void ATOMIC_ADD(value_t *target, value_t val) {
  if(val == 0) return;
  #pragma omp atomic
  *target += val;
}

inline void ATOMIC_MULT(value_t *target, value_t val) {
  if(val == 1) return;
  #pragma omp atomic
  *target *= val;
}

inline void ATOMIC_MIN(value_t *target, value_t value){
  #pragma omp critical
  {
    *target = (*target < value) ? *target : value;
  }
}

inline void ATOMIC_AND(value_t *target, value_t val) {
  #pragma omp atomic
   *target &= val;
}


inline void ATOMIC_OR(value_t *target, value_t val) {
  #pragma omp atomic
   *target |= val;
} 


#endif
