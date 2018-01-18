#ifndef NODE_INT_MAP_H
#define NODE_INT_MAP_H

#include "graph.h"

typedef struct nodeIntMapElement {
  node_t key;
  int value;
}nodeIntMapElement;

/* 
   Map data structure.
   If the map needs to be thread safe
   use nodeIntMapAtomic
*/

typedef struct nodeIntMap {
  nodeIntMapElement* list;
  int maxSize;
  int size;
} nodeIntMap;


/* 
   We uae array implementation of the map.
   Search O(n)
   Insert O(1)
   Update O(1)
 */

void __initNodeIntMap(nodeIntMap* map, int initValue) {
  map->size = 0;
  int i;
  for(i=0; i< map->maxSize; i++) {
    map->list[i].value = initValue;
  }
}

nodeIntMap* initNodeIntMap(nodeIntMap* map, int size, int initValue) {
  map = (nodeIntMap*) malloc (sizeof(nodeIntMap));
  map->maxSize = size;
  map->list = (nodeIntMapElement*) malloc (size * sizeof (nodeIntMapElement));
  __initNodeIntMap(map,initValue);
  return map;
}

nodeIntMap* reinitNodeIntMap(nodeIntMap* map, int mapSize, int initValue) {
  
  if(mapSize > map->maxSize) {
    map->maxSize = mapSize;
    map->list = (nodeIntMapElement*) realloc (map->list, mapSize * sizeof (nodeIntMapElement));    
  }
  __initNodeIntMap(map,initValue);

  return map;
}



void closeNodeIntMap(nodeIntMap* map) {
  free(map->list);
}


node_t mapMaxValueKey(nodeIntMap* map) {
  int max = 0, i;
  node_t maxVal = NIL_NODE;
  for(i=0;i<map->size;i++) {
    if(max < map->list[i].value) {
      maxVal = map->list[i].key;
      max = map->list[i].value;
    }  
  }
  return maxVal;
}


void changeValue(nodeIntMap* map, node_t key, int inc) {
  int i;
  for(i=0;i<map->size;i++) {
    if(key == map->list[i].key) {
      map->list[i].value += inc;
      break;
    }  
  }
  if(i == map->size ) {
    map->list[i].key = key;
    map->list[i].value = inc;
    map->size++;
  }
}




#ifdef ATOMICMAP
/* TODO  */
typedef struct nodeIntMapAtomic {
  nodeIntMapElement* list;
  int maxSize;
  int size;
} nodeIntMapAtomic;

void changeValueAtomicAddNodeIntMap(nodeIntMapAtomic* map, node_t key, int inc) {
#pragma omp critical
  {
    
  }
}

#endif


#endif
