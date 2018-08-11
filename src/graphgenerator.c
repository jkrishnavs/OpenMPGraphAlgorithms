/*****
 * Original base algorithm from GTGraph. 
 ****/
#include<unistd.h>
#include"graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "graphprop.h"

graph* randomGenerator();
graph* erdosRenyiGenerator();
graph* rmatGenerator();
graph* SSCAGenerator();
graph* propertyControlledGraphGenerator();
graph* callappropriategenerator();
graph* allocateMemoryforGraph(nodes_t n, edges_t m);
void doubleMergeSort(node_t* l1, node_t* l2, edge_t left, edge_t right);
void merge(node_t* l1, node_t* l2, edge_t left, edge_t mid, edge_t right);


typedef enum GraphModel {
  Random,
  ErdosRenyi,
  RMAT,
  SSCA,
  PCGG // Property Controlled Graph Generation 
} GraphModel;


typedef enum GraphProperty {
  degreeanddensity,
  clustercoeff,
  degreeaed,
  degreeanddegresd
}

GraphModel model;
boolean weighted;
nodes_t numNodes;
edges_t numEdges;
/*** Random Parameters ***/
boolean selfloop;
/***** Erdos Renyi *****/
double edgeProbability;
/***** RMAT **********/
double a,b,c,d;
/********** Graph prop value *******/
/**
 NOTE: These parameters has higher preference over numEdges and numNodes.
 if these values contradict with the numEdges and numNodes 
 we update the numNodes and numEdges based on these values.
**/
double degree;
double density;
double clusteringCoeff;
double aed;
/******
 * standard deviation on degrees
 */
double degreesd;
/*****************/
/** 
   Acceptable error as a percentage.
   the final graph properties will be
 in the range (prop*(1-err), prop*(1+err)). 
**/
double err;

// TODO device a time out mechanism

/*** Weight parameters ****/
int maxWeight;
int minWeight;




/**** dafualt configs ******/
void setbaseConfigs() {
  model = Random; // m
  weighted = false; // wf
  numNodes = 100000; //n
  numEdges = 1000000; // e
  selfloop = false; // sl
  edgeProbability = 0.0001; //ep
  a = 0.45; // a
  b = 0.15; // b
  c = 0.15; // c
  d = 0.25; // d

  /*** Properties of interest *******/
  degree  = 10;// dg
  density = 0.0001; // dn
  clusteringCoeff = 0.1; // cc
  aedbase = 0.33; // aed
  degreesd = 0.2; //dsd
  err = 1; /*1%*/ // err
  maxWeight = 100; //mw
  minWeight = 1; // iw
 
}

void updateConfigs(const char* configfile) {
  FILE *f;

  f = fopen(filename, "r");
  if(f == NULL) {
    printError(CONFIG_FILE_NOT_FOUND, 0, NULL);
    return NULL;
  }
  // we assume the first entry in the config file specifies the model
  // of graph generation.
  // TODO read configs
  
  
}



int runalgo(int argc, char** argv) {

  if(argc > 4) {
    const char* argList[3] = {" <configfile.conf>" , "<outputfile.edges>","<propfile.prop>"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }
  updateConfigs(configfile);
  graph * G = callappropriategenerator();
  assert(G != NULL);
  writeBackGraph(G, argv[2]);

  /*** Graph prop collection ***/
  avgClusterCoeff(G);
  avgEdgeDistance(G);
  diameter(G);
  sparsity(G);
  triangle_counting(G);
  printf("Avg Adjecency =  %f \n", avgincident);
  printf("\nAvg ClusterCoeff = %f \n",clusterCoeff);
  printf("Avg Edge Distance = %f \n", aed);
  printf("Sparsity = %.7f \n",sparsityMeasure);
  printf("The Scaled percentange Triangles = %f\n\n", scaledT);

  return 0;
}

graph* callappropriategenerator() {
  if(model == Random) {
    return randomGenerator();
  } else if (model  ==   ErdosRenyi) {
    return erdosRenyiGenerator();
  } else if (model == RMAT) {
    return rmatGenerator();
  } else if (model == SSCA) {
    return SSCAGenerator();
  } else if 
}



graph* randomGenerator() {
  int *stream1, *stream2;
  FILE* outfp;
    
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  graph* G = allocatememoryforGraph(numNodes, numEdges);
  G->r_node_idx = (node_t) malloc (numEdges * sizeof(node_t));

 
  
  for(edge_t i = 0; i< numEdges; i++) {
    node_t u = (node_t) isprng(stream1) % numNodes;
    node_t v = (node_t) isprng(stream1) % numNodes;
    if((u == v) && (selfloop == false)) {
      i--;
      continue;
    }
    if(weighted) {
      int w = (int)( minWeight + ((double )(maxWeight - minWeight)) * sprng(stream2));
      G->weights[i] = w;
    }
    G->node_idx[i] = v;
    G->r_node_idx[i] = u;
    G->begin[i+1]++; // array start with 0
  }
  
  for(node_t n = 1;n< numNodes; n++) {
    G->begin[n+1] += G->begin[n];
  }
  doubleMergeSort(G->r_node_idx, G->node_idx, 0, numEdges);
  return G;
}

#define comparedoubleList(l1,l2,i1,i2) ((l1[i2]< l1[i1]) ? true: ( (l1[i2]== l1[i1]) ? l2[i2] < l2[i1]: false))

#define comparedoubleList(l1,l2, r1, r2, i1, i2) ((r1[i2]< l1[i1]) ? true: ( (r1[i2]== l1[i1]) ? r2[i2] < l2[i1]: false))

inline void swapdoublelist(node_t* l1, node_t* l2, edge_t i1, edge_t i2) {
  node_t tmp;
  tmp = l1[i1]; l1[i1] = l1[i2]; l1[i2] = tmp;
  tmp = l2[i1]; l2[i1] = l2[i2]; l2[i2] = tmp;
}

inline void copydoublelist(node_t* l1, node_t* l2, node_t* c1, node_t* c2, edge_t i1, edge_t i2) {
  c1[i2] = l1[i1];
  c2[i2] = l2[i1];
}

void merge(node_t* l1, node_t* l2, edge_t start, edge_t mid, edge_t end) {
  node_t t1 = start;  node_t t2 = mid;
  node_t tmp;  node_t temp; node_t w = mid-1;
  node_t tempqueue[(mid - start)*2]; // should we use malloc ?
  node_t offset = mid-start;
  node_t *tq2 = &(tempqueue[offset]);
  
  node_t tpf = 0; node_t tpb = 0; 
  
  while(t1 < mid && t2< end) {
    if(tpf == tpb) {
      // empty queue
      if(comparedoubleList(l1, l2, t1, t2)) {
	copydoublelist(l1, l2, tempqueue, tq2, t1, tpf);
	tpf++;
	copydoublelist(l1, l2, l1, l2, t2, t1);
	t2++;
      }
    } else{
      if(comparedoublelist(tempqueue, tq2, l1, l2, tpb, t2)) {
	copydoublelist(l1, l2, tempqueue, tq2, t1, tpf);
	tpf++;
	copydoublelist(l1, l2, l1, l2, t2, t1);
	t2++;
      } else {	
	copydoublelist(l1, l2, tempqueue, tq2, t1, tpf);
	tpf++;
	copydoublelist(tempqueue, tq2, l1, l2,  tpb, t1);
	tpb++;
      }
    }
    t1++;
  }
  if(t1 < mid) {
    // on highly rare occations
    // copy rest of the first half to the temp array
    while(t1 < mid) {
      copydoublelist(l1, l2, tempqueue, tq2, t1, tpf);
      tpf++;	
    }
    // now copy back withou comparison t1 array already sorted.
    while(tpf > tpb ) {
      copydoublelist(tempqueue, tq2,  l1, l2,  tpf, t1);
      t1++; tpb++;
    }
  } else {
    while(tpf > tpb ) {
      if(t2 < end && comparedoublelist(tempqueue, tq2, l1, l2, tpb, t2)) {
	copydoublelist(l1, l2, l1, l2, t2, t1);
	t2++;
      } else {
	copydoublelist(tempqueue, tq2, l1, l2,  tpb, t1);
	tpb++;
      }
      t1++;
    }
  }
  // TODO add assert?
  assert(tpf == tpb);
  assert(t1 == (t2-1)); 
}
 
 



void doubleMergeSort(node_t* l1, node_t* l2, edge_t left, edge_t right) {

  if((right - left)< 100) {
    edge_t mid = left + (right-left)/2;
    doubleMergeSort(l1,l2, left, mid);
    doubleMergeSort(l1,l2, mid, right);
    merge( l1, l2, left, mid, right);
  } else {
    edge_t itr = left;
    for(; itr < right; itr++) {
      edge_t itr2= itr;
      for(; itr2 < right; itr2++) {
	if(comparedoubleList(l1,l2, itr, itr2)) {
	  swapdoubleList(l1,l2, itr, itr2);
	}
      }
    }
  }

}



graph* allocateMemoryforGraph(nodes_t n, edges_t m) {
  graph *G = (graph*) malloc (sizeof(graph));
  G->begin = (edge_t*) malloc( (n+1)* sizeof(edge_t));
  G->node_idx = (node_t*) malloc(m * sizeof(node_t));
  if(weighted == true) {
    G->weights = (int*) malloc(m* sizeof(int));
  }
  /* Initialize */
  memset(G->begin, 0, (n+1) * sizeof(edge_t));
  
  return G;
}


graph* erdosRenyiGenerator() {
  int *stream1, *stream2;
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  edge_t maxEdges = 0;
  numEdges = 0;
  if(selfloop == false) {
    /* Should cover 90 percent of edges */
    maxEdges = (edge_t) ( 1.1 * ((edgeProbablity) * numNodes) * (numNodes -1));
  } else {
    maxEdges = (edge_t) ( 1.1 * ((edgeProbablity) * numNodes) * (numNodes));
  }
  graph *G = allocatememoryforGraph(n,maxEdges);
   
  /* Write the no. of edges later */

  edg = 0;
  for (node_t i=0; i<n; i++) {
    G->begin[i] = edg;
    for (node_t j=0; j<n; j++) {
      if ((i==j) && (selfloop == false))		
	continue;
      if (p > sprng(stream1)) {
	if(weighted) {
	  int w = (int)( minWeight + ((double )(maxWeight - minWeight)) * sprng(stream2));
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  assert(edg < maxEdges);
  return G;
}

graph* rmatGenerator() {
  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/

}

graph* propertyControlledGraphGenerator() {


}


graph* SSCAGenerator() {
    /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/

  // TODO
}
