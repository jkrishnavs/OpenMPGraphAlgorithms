#include<stdlib.h>
#include "graph.h"
#include "mainFunctions.h"
#include "print.h"
//#include "powerperformacetracking.h"
//#include "communities.h"
#include "graphprop.h"
#include "nodeIntMap.h"


#define DEBUG_ON



int adj = 0; // 1 for reverse adjlist or 0 in adjlist.
bool skewed  =  false;

node_t* comm = NULL;

typedef struct graphmap {
  node_t comm;
  node_t newPos;
  node_t revPos;
} graphmap;


typedef struct ClusterDetails {
  int id;
  int numNodes;
  int numEdges;
  int external;
  int commDistance;
  //  node_t* nodeList;
} ClusterDetails;

inline void CopyClusterDetails(ClusterDetails* cd,int src, int dest) {
  cd[dest].id = cd[src].id;
  cd[dest].numNodes = cd[src].numNodes;
  cd[dest].numEdges = cd[src].numEdges;
  cd[dest].external = cd[src].external;
  cd[dest].commDistance = cd[src].commDistance;  
}



graph* createPreprocessedGraph(graph *G, graphmap * gm) {
  /* Create the copy */
  graph* newG = createGraph();

  bool hasEdgeWeight= false;
  

  newG->numNodes = G->numNodes;
  newG->numEdges = G->numEdges;
  newG->begin = (edge_t*) malloc (sizeof(edge_t) * (newG->numNodes+1));
#ifdef DEBUG_ON
  assert(newG->begin != NULL);
#endif
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
commlist = nodes in community
stacklist =  stack 
visited = visited list
commid community id
 */
void * topologicalSort(graph * G, node_t* commlist , node_t* stacklist, int* visited,  node_t* comm, node_t* commpos, int commSize, node_t commId) {
    int sortedSize = 0;
    stacklist[sortedSize] = commlist[0];
    visited[0] = 1;
    sortedSize++;
    int curIndex = 0;
    node_t source;
    // printf("The start of stack \n");
    while(sortedSize < commSize) {
      if(curIndex == sortedSize) {
	for(int i =1; i< commSize; i++ ) {
	  if(visited[commpos[commlist[i]]] == 0) {
	    stacklist[sortedSize] = commlist[i];
	    visited[commpos[commlist[i]]] = 1;
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
	if(comm[d] == commId && visited[commpos[d]] == 0) {
	  stacklist[sortedSize] = d;
	  visited[commpos[d]] = 1;
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

void initCommunities(graph *G) {
  
  comm = (node_t*) malloc (G->numNodes * sizeof(node_t));
  assert(comm != NULL);

}

double* coeff;

#define cohval(index,arr) arr[index]

void merge_serial(node_t* index,  double* vals, node_t start, node_t mid, node_t end) {
  node_t t1 = start;  node_t t2 = mid;
  node_t tmp;  node_t temp; node_t w = mid-1;
  node_t tempqueue[mid - start]; // should we use malloc ?
  node_t tpf = 0; node_t tpb = 0; 
  
  while(t1 < mid && t2< end) {
    if(tpf == tpb) {
      // empty queue
      if(cohval(index[t1], vals) < cohval(index[t2], vals)) {
	tempqueue[tpf] = index[t1];
	tpf++;
	index[t1] = index[t2];
	t2++;
      }
    } else{
      if(cohval(tempqueue[tpb], vals) < cohval(index[t2], vals)) {
	tempqueue[tpf] = index[t1];
	tpf++;
	index[t1] = index[t2];
	t2++;
      } else {
	tempqueue[tpf] = index[t1];
	tpf++;
	index[t1] = tempqueue[tpb];
	tpb++;
      }
    }
    t1++;
  }
  if(t1 < mid) {
    // on highly rare occations
    // copy rest of the first half to the temp array
    while(t1 < mid) {
      tempqueue[tpf] = index[t1];
      tpf++;	
    }
    // now copy back withou comparison t1 array already sorted.
    while(tpf > tpb ) {
      index[t1] = tempqueue[tpb];
      t1++; tpb++;
    }
  } else {
    while(tpf > tpb ) {
      if(t2 < end && cohval(tempqueue[tpb], vals) < cohval(index[t2], vals)) {
	index[t1] = index[t2];
	t2++;
      } else {
	index[t1] = tempqueue[tpb];
	tpb++;
      }
      t1++;
    }
  }
  // TODO add assert?
  assert(tpf == tpb);
  if(t2 < end) {
    assert(t1 == (t2-1));
  } else {
    //assert(t2 == end);
    assert(t1 == end);
  }
  //
  
}
void merge_parallel(node_t* index1,  double* vals, node_t start, node_t mid, node_t end) {
  // TODO
}
void sort_selection(node_t* index, double* vals, node_t start, node_t end) {
  node_t sindex;
  for(node_t i = start; i< end-1; i++) {
    sindex = i;
    double coh = cohval(index[i], vals);
    for(node_t j = start+1; j<end; j++) {
      double coj = cohval(index[j], vals);
      if(coj > coh) {
	sindex = j;
	coh = coj;
      }
    }
    node_t temp = index[i];
    index[i] = index[sindex];
    index[sindex] = temp;
  }
}


void sort_serial(node_t* index, double* vals, node_t start, node_t end) {
  if((end-start) < 1024) {
    sort_selection(index, vals, start, end);
  } else {
    node_t midpoint = (start+end)/2;
    sort_serial(index, vals, start, midpoint);
    sort_serial(index, vals, midpoint, end);
    merge_serial(index, vals, start, midpoint, end);
  }
}

void sort_parallel(node_t* index, double* vals, node_t start, node_t end) {

  if(end - start < 8192) {
    sort_serial(index, vals, start, end);
  }
  else {
/* #pragma omp parallel */
/*     { */
#pragma omp task 
      sort_parallel(index, vals, start , (start+end)/2);
#pragma omp task 
      sort_parallel(index, vals, (start+end)/2 , end);
#pragma omp taskwait
      merge_serial(index, vals, start, (start+end)/2, end);
    /* } */
  }
}

void sort(node_t* index, double* vals, node_t start, node_t end) {
  if((end - start) < 8192) {
    sort_serial(index, vals,start,end);
  } else {
#pragma omp parallel
    {
      sort_parallel(index, vals, start, end);
    }
  }
}

/* int comp (const void * elem1, const void * elem2) { */
/*   double f = coeff[*((int*)elem1)]; */
/*   double s = coeff[*((int*)elem2)]; */
/*   if (f > s) return -1; */
/*   if (f < s) return 1; */
/*   return 0;   */
/* } */





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
	  //if((G->begin[x0+1] - G->begin[x0]) > 1) {
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
	    //}
	  
	}
	closeNodeIntMap(map);
      }
      itrs++;
    } while ( !finished && maxItrs > itrs);

  gettimeofday(&end, NULL);

  
  printf("Community detection The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  gettimeofday(&start, NULL);
  
  ClusterDetails* cd = (ClusterDetails*) malloc (G->numNodes * sizeof(ClusterDetails));
  graphmap* gm = (graphmap*) malloc (G->numNodes * sizeof(graphmap));


#pragma omp parallel for
  for(node_t i=0; i<G->numNodes; i++ ) {
    gm[i].newPos = NIL_NODE;
    gm[i].revPos = NIL_NODE;
    cd[i].numNodes = 0;
    cd[i].numEdges = 0;
    cd[i].external = 0;
    cd[i].id = i;
    //cd[i].nodeList = NULL;
  }


  /* Position of the node in the community */
  node_t *commpos = (node_t*) malloc (G->numNodes * sizeof(node_t));
  
#ifdef DEBUG_ON
  printf("Start the data collection process \n");
  /**
   * IMPORTANT: There is a high probability that in the input graph
   * the adject nodes are in the same cluster.
   * As a result, we first compare the previous nodes's 
   * cluster ID.
   **/
  
#endif

  for(node_t i=0; i< G->numNodes; i++) {
    node_t comid = comm[i];
    commpos[i] = cd[comid].numNodes;
    cd[comid].numNodes++;
    for(edge_t e = G->begin[i]; e < G->begin[i+1]; e++) {
      node_t end = G->node_idx[e];
      if(comm[end] == comid)
	cd[comid].numEdges++;
      else
	cd[comid].external++;
    }
  }
  
/* number of communities */
 node_t noofComm = 0;
 for(node_t i=0; i< G->numNodes; i++) {
   if(cd[i].numNodes > 0) {
     if(noofComm < i) {
       CopyClusterDetails(cd, i, noofComm);
     }
#ifdef DEBUG_ON
     assert(cd[i].numNodes == cd[noofComm].numNodes);
     //     printf("The node id %d numNodes %d  i val is %d and i numNodes is %d \n", noofComm, cd[noofComm].numNodes, i , cd[i].numNodes);
     assert(cd[i].numNodes > 0);
     // exit(0);
#endif
     noofComm++;     
   }
 }
  

#ifdef DEBUG_ON
  printf("The number of clusters is %d \n ", noofComm);
  /* for(node_t debclusters = 0; debclusters < noofComm; debclusters++) { */
  /*   printf("The cluster id is %d and number of nodes in cluster is %d \n", cd[debclusters].id, cd[debclusters].numNodes); */
  /* } */
#endif

  

  gettimeofday(&end, NULL);
  printf("data collection The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
    

  /* 
     The community with highest diffrenece between internal nodes 
     and external nodes per node 
     will be psuhed to start.
    */
  gettimeofday(&start, NULL);

  /*community positions*/
  int* commDetPointers  = (int*) malloc (noofComm * sizeof(int));
  coeff = (double*) malloc (noofComm * sizeof(double));
  

#pragma omp parallel for
  for(node_t i=0; i < noofComm; i++) {
#ifdef DEBUG_ON
    // printf("The node id %d numNodes %d \n",i , cd[i].numNodes);
    assert(cd[i].numNodes > 0);
#endif
    commDetPointers[i] = i;
    coeff[i] = ((double)cd[i].numEdges)/cd[i].numNodes;
  }
  
  sort_serial(commDetPointers, coeff, 0, noofComm);
  free(coeff);
  gettimeofday(&end, NULL);
  printf("sorting The required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));
  /* Decide on destinations  */
  gettimeofday(&start, NULL);
  node_t commEnd = 0;
  node_t commStart = 0; 	/* These are the start and end indexes of graph Map */  
  for(node_t c =0; c<noofComm; c++) {
    node_t cid = commDetPointers[c];
    cd[cid].commDistance = commStart;
    assert(cd[cid].numNodes > 0);
    commStart += cd[cid].numNodes;
  }
  assert(commStart == G->numNodes);

  /** Collect the nodes of a cluster together **/
  node_t* clusterednodes = (node_t*) malloc (G->numNodes * sizeof(node_t));
  node_t* clusterPositions = (node_t*) malloc (G->numNodes * sizeof(node_t));
#pragma omp parallel for
  for(node_t c=0; c< noofComm;c++) {
    clusterPositions[cd[c].id] = c;
  }
#ifdef DEBUG_ON
#pragma omp parallel for
  for(node_t i =0;i< G->numNodes;i++) {
    clusterednodes[i] = NIL_NODE;
  }
  for(node_t c; c< noofComm; c++) {
    assert(clusterPositions[cd[c].id] == c); 
  }
#endif  
  for(node_t i =0;i< G->numNodes;i++) {
    node_t nodeCluster = comm[i];
    node_t clusterPos = clusterPositions[nodeCluster];
    node_t nodePos = cd[clusterPos].commDistance + commpos[i];
#ifdef DEBUG_ON
    //    printf("The node id %d community id %d (%d), community start is %d nodePos =  %d clusterednode %d \n",i,comm[i], cd[clusterPos].id, cd[clusterPos].commDistance, nodePos, clusterednodes[nodePos]);
    assert(cd[clusterPos].id == comm[i]);
    assert(clusterednodes[nodePos] == NIL_NODE);
#endif
    clusterednodes[nodePos] = i;
  }

  free(clusterPositions);
  
  
#pragma omp parallel
  {
    int* visitedlist = (int*) malloc (cd[0].numNodes * sizeof(int));
    node_t * commlist = (node_t*) malloc (cd[0].numNodes * sizeof(node_t));  
    node_t* stacklist = (node_t*) malloc (cd[0].numNodes * sizeof(node_t));
    node_t commMax = cd[0].numNodes;

#pragma omp for
    for(node_t c = 0; c< noofComm; c++) {
      node_t cid = commDetPointers[c];
      node_t myCommunitySize = cd[cid].numNodes;
      node_t commId = cd[cid].id;
      node_t myStart = cd[cid].commDistance;
      if(myCommunitySize > commMax) {
	commlist = (node_t*) realloc(commlist, myCommunitySize * sizeof(node_t));
	stacklist = (node_t*) realloc(stacklist, myCommunitySize * sizeof(node_t));
	visitedlist = (int*) realloc(visitedlist , myCommunitySize * sizeof(int));
	commMax = myCommunitySize;
      }
      for(node_t i = 0; i< myCommunitySize; i++) {
	commlist[i] = clusterednodes[i+ myStart];
#ifdef DEBUG_ON
	assert(comm[commlist[i]]  == commId);
#endif
      }
      /*************************** TODO *************************************/
      if(myCommunitySize > 3) {
	memset(visitedlist, 0,myCommunitySize  * sizeof(int));
	// BFS
	topologicalSort(G, commlist, stacklist, visitedlist, comm, commpos, myCommunitySize, commId);
	/* TODO	   Do a BFS after reversing the list.	*/
	/* We might not require this exchange after the TODO */
	node_t* temp = commlist;
	commlist = stacklist;
	stacklist = temp;
      }
      /*  The order is set in commlist */
      for(node_t i =0;i < myCommunitySize; i++){
	assert(gm[myStart +  i].revPos == NIL_NODE);
	assert(gm[commlist[i]].newPos == NIL_NODE);
	gm[myStart +  i].revPos = commlist[i];
	gm[commlist[i]].newPos = myStart + i;
      }
    }
    free(stacklist);
    free(commlist);
    free(visitedlist);
  }
  free(comm);
  free(commpos);
  free(cd);
  free(commDetPointers);
  free(clusterednodes);

  gettimeofday(&end, NULL);
  printf("Internal sorting the required time is %f \n",((end.tv_sec - start.tv_sec)*1000 + ((double)(end.tv_usec - start.tv_usec))/1000));

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
  if(argc < 7) {
    const char* argList[6] = {" <inputfile> " , "graphformat.txt","<outputfile>", "<outputpropfile>", "<outputmapfilename>", "<adjecencyflag>"};
    printError(INCORRECT_ARG_LIST, 6, argList);
    return -1;

  }
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










