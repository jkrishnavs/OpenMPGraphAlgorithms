#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "communities.h"
#include "graphprop.h"



int adj = 0; // 1 for reverse adjlist or 0 in adjlist.
bool skewed  =  false;


typedef struct graphmap {
  node_t comm;
  node_t newPos;
  node_t revPos;
} graphmap;


typedef struct CommunityDetails {
  int id;
  int numNodes;
  int numEdges;
  int external;
  int commDistance;
} CommunityDetails;


graph* createPreprocessedGraph(graph *G, graphmap * gm) {
  /* Create the copy */
  graph* newG = createGraph();

  bool hasEdgeWeight= false;
  

  newG->numNodes = G->numNodes;
  newG->numEdges = G->numEdges;
  newG->begin = (edge_t*) malloc (sizeof(edge_t) * (newG->numNodes+1));
  assert(newG->begin != NULL);
  newG->node_idx = (node_t*) malloc(sizeof(node_t) * newG->numEdges);
  edge_t edgepos = 0;
  newG->begin[0] = 0;
  for(int i=0;i< newG->numNodes; i++) {
    assert(gm[i].revPos != NIL_NODE);
    node_t origPos = gm[i].revPos;
    // reserse the edge list size
    newG->begin[i+1] = newG->begin[i] + (G->begin[origPos+1] - G->begin[origPos]);
    edge_t st = newG->begin[i];
    edge_t ed = newG->begin[i];

    // add edges
    for(edge_t e = G->begin[origPos]; e < G->begin[origPos+1]; e++) {
      node_t end = G->node_idx[e];
      assert(gm[end].newPos != NIL_NODE);
      // get from map the new node id
      node_t newEnd = gm[end].newPos;
      edge_t sorte = st;
      // find right position
      for(; sorte < ed ; sorte++) {
	if(newEnd < newG->node_idx[sorte])
	  break;
      }
      // shift edges
      for(edge_t sf = ed-1; sf >= sorte; sf --) {
	newG->node_idx[sf+1] = newG->node_idx[sf];  
      }
      // add to right position
      newG->node_idx[sorte] = newEnd;
      ed++;
    }
    assert(ed == newG->begin[i+1]);
  }
  return newG;
}



/* 
   TODO there can be atmost 1 node with no incomming edges in the community. 
   (The universal source node of the community).
   If we can somehow find that node (if present, possibly it will the 
   community ID TODO VERIFY AND OPTIMIZE THIS), the we can get away with creating the stacklist
   and commlist and directly use the G and comm to create the consolidated list.
 */
void * partialBFS(graph * G, node_t* commlist , node_t* stacklist, int* visited,  node_t* comm, int commSize, node_t commId) {
    int sortedSize = 0;
    stacklist[sortedSize] = commlist[0];
    visited[commlist[0]] = 1;
    sortedSize++;
    int curIndex = 0;
    node_t source;
    // printf("The start of stack \n");
    while(sortedSize < commSize) {
      if(curIndex == sortedSize) {
	for(int i =1; i< commSize; i++ ) {
	  if(visited[commlist[i]] == 0) {
	    stacklist[sortedSize] = commlist[i];
	    visited[commlist[i]] = 1;
      	    sortedSize++;
	    break;
	  }
	}
      }
      assert(curIndex< sortedSize);
      source = stacklist[curIndex];
      curIndex++;
      for(edge_t e = G->begin[source]; e < G->begin[source+1]; e++) {
	node_t d = G->node_idx[e];
	if(comm[d] == commId && visited[d] == 0) {
	  stacklist[sortedSize] = d;
	  visited[d] = 1;
	  sortedSize++;
	}
      }
    }
    //    printf("The end of stack \n");
    assert(sortedSize == commSize);
}


void dumpmapping(graph *G, graphmap*  gm,const char * filename) {
  FILE *fp = fopen(filename, "w");
  for(node_t i = 0; i < G->numNodes; i++) {
    fprintf(fp, "%d\t%d\n", i, gm[i].newPos);
  }
}


graph* preprocess(graph* G, const char* mapfile) {
  bool finished = false ;
  finished = true ;
  double mean = ((double)G->numEdges)/ G->numNodes;
  double upperlimit = 1.5 * mean;
  double lowerlimit = mean * 0.5;
  // int outliers = 0;


  bool hasEdgeWeight = false;

  if(G->weights!= NULL)
    hasEdgeWeight = true;

  
  struct timeval start, end;
  
  if(adj == 1) {
    /* reverse edge parallelism */
    createReverseEdges(G);
    edge_t* r_begin = G->r_begin;
    node_t* r_node_idx = G->r_node_idx;
    free(G->begin);
    free(G->node_idx);
    G->begin = r_begin;
    G->node_idx = r_node_idx;
  }


  gettimeofday(&start, NULL);
  

  initCommunities(G);

#pragma omp parallel
  {
#pragma omp for
    for (node_t x = 0; x < G->numNodes; x ++)
      comm[x] = x ;
  }

  int maxItrs = 100;
  
  int itrs = 0;
  do
    {
      finished = true ;
#pragma omp parallel
      {
	nodeIntMap *map;
	map = NULL;
	map = initNodeIntMap(map, 32, 0);
	node_t x0;
#pragma omp for schedule(dynamic, PAR_CHUNKSIZE)
	for (x0 = 0; x0 < G->numNodes; x0 ++) {
	  /* We classify only nodes with more than one edge.
	     The non classified nodes will be pused to last. */
	  if((G->begin[x0+1] - G->begin[x0]) > 1) {
	    map = reinitNodeIntMap(map, G->begin[x0+1] - G->begin[x0], 0);
	    for (edge_t y_idx = G->begin[x0];y_idx < G->begin[x0+1] ; y_idx ++) {
	      node_t y = G->node_idx [y_idx];
	      node_t source;
	      source = comm[y];
	      changeValue(map, source, 1);
	    }
	    node_t maxVal = mapMaxValueKey(map);
	    if ( comm[x0] != maxVal  && maxVal != NIL_NODE) {
#pragma omp atomic write
	      comm[x0] = maxVal;
	      finished = false ;
	    }
	  }
	  
	}
	closeNodeIntMap(map);
      }
      itrs++;
    } while ( !finished && maxItrs > itrs);

  gettimeofday(&end, NULL);

  
  printf("Community detection The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  
  CommunityDetails* cd = (CommunityDetails*) malloc (G->numNodes * sizeof(CommunityDetails));
  graphmap* gm = (graphmap*) malloc (G->numNodes * sizeof(graphmap));
#pragma omp parallel for
  for(node_t i=0; i<G->numNodes; i++ ) {
    gm[i].newPos = NIL_NODE;
    gm[i].revPos = NIL_NODE;
    cd[i].numNodes = 0;
    cd[i].numEdges = 0;
    cd[i].external = 0;
    cd[i].id = (int)NIL_NODE;
  }


  /* Number of communities */
  node_t noofComm = 0;
  gettimeofday(&start, NULL);
  for(node_t i=0; i< G->numNodes; i++) {
    node_t comid = comm[i];
    node_t index = NIL_NODE;
    for(node_t j =0; j <noofComm; j++) {
      if(cd[j].id == comid) {
	index = j;
	cd[index].numNodes++;
	break;
      }
    }
    if(index == NIL_NODE) {
      index = noofComm; 
      cd[index].id = comid;
      cd[index].numNodes = 1;
      noofComm++;
    }
    for(edge_t e = G->begin[i]; e < G->begin[i+1]; e++) {
      node_t end = G->node_idx[e];
      if(comm[end] == comid)
	cd[index].numEdges++;
      else
	cd[index].external++;
    }
  }


  gettimeofday(&end, NULL);
  printf("data collection The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
    

  /* The community with highest diffrenece between internal nodes 
     and external nodes per node 
    will be psuhed to start.
    */

  gettimeofday(&start, NULL);
  
  for(node_t i=0; i< noofComm -1; i++) {
    for(node_t j= i+1; j < noofComm; j++) {
      double coh = ((double)cd[i].numEdges)/cd[i].numNodes;
      double coj = ((double)cd[j].numEdges)/cd[j].numNodes;
      if(coj > coh) {
	/* exchage the nodes */
	CommunityDetails temp;
	memcpy(&temp, &cd[i], sizeof(CommunityDetails));
	memcpy(&cd[i], &cd[j], sizeof(CommunityDetails));
	memcpy(&cd[j], &temp, sizeof(CommunityDetails));	
      }
    }
  }

  gettimeofday(&end, NULL);
  printf("sorting The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  /* Decide on destinations  */

  node_t * commlist = (node_t*) malloc (cd[0].numNodes * sizeof(node_t));
  node_t* stacklist = (node_t*) malloc (cd[0].numNodes * sizeof(node_t));

  
  int * visitedlist = (int*) malloc (G->numNodes * sizeof(int));

  node_t commStart = 0; 	/* These are the start and end indexes of graph Map */
  node_t commEnd = 0;

  node_t commMax = cd[0].numNodes;


  gettimeofday(&start, NULL);

  for(node_t c = 0;c < noofComm; c++) {
    node_t commSize = cd[c].numNodes;
    //    printf("The community Size is %d internal edges %d \n", commSize, cd[c].numEdges);
    node_t commId = cd[c].id;
    if(commSize > commMax) {
      commlist = (node_t*) realloc(commlist, commSize * sizeof(node_t));
      stacklist = (node_t*) realloc(stacklist, commSize * sizeof(node_t));

      commMax = commSize;
    }
    
    node_t t = 0;
    /* Get all nodes in the community */
    for(node_t i =0; i< G->numNodes; i++) {
      if(comm[i] == commId) {
	commlist[t] = i;
	t++;
      }
    }

    if(t != commSize) {
      printf("The values is %d %d \n", t, commSize);
      // set visited to 0;
    }
    assert(t == commSize);


    if(commSize > 3) {
      memset(visitedlist, 0, G->numNodes *sizeof(int));
      // BFS
      partialBFS(G, commlist, stacklist, visitedlist, comm, commSize, commId);

      /* TODO
	 Do a BFS after reversing the list.
      */

      /* We might not require this exchange after the TODO */
      node_t* temp = commlist;
      commlist = stacklist;
      stacklist = temp;
	
      /* memset(visitedlist, 0, G->numNodes *sizeof(int)); */
      /* partialBFS(G, stacklist, commlist, visitedlist, commSize, commId); */
    }
      /*  The order is set in commlist */
    for(node_t i =0;i < commSize; i++){
      assert(gm[commStart +  i].revPos == NIL_NODE);
      assert(gm[commlist[i]].newPos == NIL_NODE);
      gm[commStart +  i].revPos = commlist[i];
      gm[commlist[i]].newPos = commStart + i;
    }
    commStart += commSize; 
  }
  
  
  gettimeofday(&end, NULL);
  printf("Internal sorting the required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));

  assert(commStart == G->numNodes);

  free(stacklist);
  free(commlist);
  free(visitedlist);


  gettimeofday(&start, NULL);
  graph* newG = createPreprocessedGraph(G, gm);
  gettimeofday(&end, NULL);
  printf("Copy back the required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  
  
  if(adj == 1) {
    /* Get the reverse edges */
    createReverseEdges(newG);
    edge_t* r_begin = newG->r_begin;
    node_t* r_node_idx = newG->r_node_idx;
    free(newG->begin);
    free(newG->node_idx);
    newG->begin = r_begin;
    newG->node_idx = r_node_idx;
  }
  
  free(cd);
  free(comm);
  dumpmapping(G, gm, mapfile);


  /*now update edgeWeights in newG*/
  if(G->weights != NULL) {
    int w = 0;
    bool found  = false;
    edge_t pos;
    newG->weights = (int*) malloc(sizeof(int) * newG->numEdges);
    node_t x0;
    for (x0 = 0; x0 < G->numNodes; x0 ++) {
      for (edge_t y = G->begin[x0];y < G->begin[x0+1] ; y ++) {
	w = G->weights[y];
	node_t d =  G->node_idx[y];      
	node_t newS = gm[x0].newPos;
	node_t newD = gm[d].newPos;
	/* assert edge is added */
	found = false;
	pos = NIL_EDGE;
	for(edge_t newY = newG->begin[newS]; newY < newG->begin[newS+1]; newY++) {
	  if(newG->node_idx[newY] == newD) {
	    found = true;
	    pos = newY;
	    break;
	  }
	}
	assert(pos != NIL_EDGE);
	newG->weights[pos] = w;
      }
    }
  }
  
  free(gm);
  
  
  return newG;
}

double avgincident;
double scaledT;
void writeSchema(const char *filename) {
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "Avg Adjecency =  %f \n", avgincident);
  fprintf(fp, "Avg ClusterCoeff = %f \n",clusterCoeff);
  fprintf(fp, "Avg Edge Distance = %f \n", aed);
  fprintf(fp, "Sparsity = %0.7f \n",sparsityMeasure);
  fprintf(fp,"The Scaled percentange Triangles = %f\n\n", scaledT);  
  fclose(fp);
}

/***
 * Common entry point for all algorithms,
 **/
int runalgo(int argc,char** argv) {
  graph* G = readGraph(argv[1], argv[2]);
  adj = atoi(argv[6]);
  graph* newG = preprocess(G, argv[5]);

  avgincident = ((double) G->numEdges )/ G->numNodes;
  
  writeBackGraph(newG, argv[3]);
  avgClusterCoeff(newG);
  avgEdgeDistance(newG);
  diameter(newG);
  sparsity(newG);
  triangle_counting(newG);
  scaledT = ((double) T )/ G->numNodes; 
  sccIndex(newG); 
  writeSchema(argv[4]);
  return 0;
}

inline void kernel(graph *G) {

}










