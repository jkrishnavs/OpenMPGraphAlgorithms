inline void ATOMIC_ADD(value_t &target, value_t val) {
  if(value == 0) return;
  #pragma omp atomic
  target += value;
}

inline void ATOMIC_INC(value_t &target) {
  #pragma omp atomic
  target ++;
}

inline void ATOMIC_MULT(value_t &target, value_t val) {
  if(value == 1) return;
  #pragma omp atomic
  target *= val;
}

inline void ATOMIC_MIN(value_t &target, value_t value){
  #pragma omp critical
  {
    target = (target < value) ? target : value;
  }
}

inline void ATOMIC_AND(value &target, value_t val) {
  #pragma omp atomic
  target &= val;
}


inline void ATOMIC_OR(value &target, value_t val) {
  #pragma omp atomic
  target |= val;
} 
