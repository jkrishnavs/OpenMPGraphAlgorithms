/*****
 * Original base algorithm from GTGraph. 
 ****/

#include<unistd.h>
#include"graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "graphprop.h"
#include "normaldistribution.h"


// TODO double check random number generation.




graph* allocateMemoryforGraph(nodes_t n, edges_t m);
void doubleMergeSort(node_t* l1, node_t* l2, edge_t left, edge_t right);
void merge(node_t* l1, node_t* l2, edge_t left, edge_t mid, edge_t right);
void choosePartition(node_t* u, node_t* v, node_t step, int* stream1);
void varyParams(double* a, double* b, double* c, double* d, int * stream3, int* stream4);


typedef enum GraphModel {
  Random,
  ErdosRenyi,
  RMAT,
  SBM, // Stochastic Block Model
  SSBM, // Symmetric Stochastic Block Model
  LSSBM, // Log Sum Stochastic Block Model
  PCGG // Property Controlled Graph Generation 
} GraphModel;



// TODO double or float for individual degree, coeff and ed ?
typedef struct MetaData {
  int changes;
  int negativeChanges;
  node_t maxDegree;
  double degree;
  double aed;
  double density;
  double degreesd;
  double clustercoeff; 
  //  node_t* degree;
  double* coeff;
  //  node_t* ed;
}gdata;


typedef enum GraphPropertyFlag {
  degreeanddensity,
  clustercoeff,
  degreeaed,
  degreeanddegresd
}GraphPropertyFlag;


typedef struct GraphProperty{
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
  /******** SSCA ********/
  node_t maxCliqueSize;
  node_t totClusters;
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
  /**
   ** If set to true the nodes will be 
   **/
  bool sortedSBMGraph; 
  
  bool degreeflag;
  bool densityflag;
  bool clusteringCoeffflag;
  bool aedflag;
  bool degreesdflag;
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
} GraphProperty;




/**** dafualt configs ******/
void setbaseConfigs(GraphProperty& p) {
  p.weighted = false; // wf
  p.selfloop = false; // sl
  p.sortedSBMGraph = true; 

  p.degreeflag = false;
  p.densityflag = false;
  p.clusteringCoeffflag = false;
  p.aedflag = false;
  p.degreesdflag = false;


  p.model = Random; // m
  p.numNodes = 100000; //n
  p.numEdges = 1000000; // e
  p.edgeProbability = 0.0001; //ep
  
  p.a = 0.45; // a
  p.b = 0.15; // b
  p.c = 0.15; // c
  p.d = 0.25; // d


  /**** Ma Clique Size **********/
  p.maxCliqueSize = 4;

  /*** Properties of interest *******/
  p.degree  = 10;// dg
  p.density = 0.0001; // dn
  p.clusteringCoeff = 0.1; // cc
  p.aedbase = 0.33; // aed
  p.degreesd = 0.2; //dsd
  p.err = 1; /*1%*/ // err
  p.maxWeight = 100; //mw
  p.minWeight = 1; // iw
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

graph* randomGenerator(const GraphProperty& p);
graph* erdosRenyiGenerator(const GraphProperty& p);
graph* rmatGenerator(const GraphProperty& p);
graph* propertyControlledGraphGenerator(GraphProperty& p);

graph* callappropriategenerator(GraphPorperty& p);
graph* StocasticBlockModel(GraphProperty& p);
graph* lowdegreeNetworkgraph(const GraphProperty& p);
graph* diagonalGraphGenerator(boolean randomizedmidpoint, const GraphPorperty& p);





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



graph* randomGenerator(const GraphProperty& p) {
  int *stream1, *stream2;
  FILE* outfp;
    
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  graph* G = allocatememoryforGraph(p.numNodes, p.numEdges);
  G->r_node_idx = (node_t) malloc (numEdges * sizeof(node_t));

 
  
  for(edge_t i = 0; i< p.numEdges; i++) {
    node_t u = (node_t) isprng(stream1) % p.numNodes;
    node_t v = (node_t) isprng(stream1) % p.numNodes;
    if((u == v) && (p.selfloop == false)) {
      i--;
      continue;
    }
    if(weighted) {
      int w = (int)( p.minWeight + ((double )(p.maxWeight - p.minWeight)) * sprng(stream2));
      G->weights[i] = w;
    }
    G->node_idx[i] = v;
    G->r_node_idx[i] = u;
    G->begin[i+1]++; // array start with 0
  }
  
  for(node_t n = 1;n< p.numNodes; n++) {
    G->begin[n+1] += G->begin[n];
  }
  doubleMergeSort(G->r_node_idx, G->node_idx, 0, p.numEdges);
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


/*** Update Graph Meta data **/
void updateGdata(graph* G, gdata& data) {
  data.changes = 0;
  data.negativeChanges = 0;
  data.maxDegree = 0;

  int maxVal = 0;
  long long degreesum = 0;

#pragma omp parallel for reduction(+: degreeSum) reduction(max: maxVal)
  for (v = 0; v < G->numNodes; v ++) {
    if((G->begin[v+1] - G->begin[v]) > maVal) {
      maxVal = (G->begin[v+1] - G->begin[v]);
    }
    degreeSum += (G->begin[v+1] - G->begin[v]);
  }

  data.maxDegree = maxVal;
  data.degree = (double) degreeSum/G->numNodes;
  data.density = (double) degreeSum/ ((G->numNodes)* (G->numNodes -1));
  double degreesd = 0;
  double xval;
    
#pragma omp parallel for reduction(+: degreesd) private(xval)
  for (v = 0; v < G->numNodes; v ++) {
    xval = (G->begin[v+1] - G->begin[v])
      degreesd +=  (xval -data.degree) * (xval -data.degree) ;
  }
  data.degreesd = sqrt(degreesd/ G->numNodes);
  data.aed = avgEdgeDistance(G);
	
  
  data.coeff = (double*) malloc (sizeof(double) * G->numNodes);
  /** coeff **/  
#pragma omp parallel
  {
    node_t v;
#pragma omp  for schedule(static)
    for (v = 0; v < G->numNodes; v ++) {
      edge_t u_idx;
      data.coeff[v] = 0;
      for (u_idx = G->begin[v]; u_idx < G->begin[v+1]; u_idx ++) {
	node_t u = G->node_idx [u_idx];
	edge_t w_idx;
	for (w_idx = u_idx+1; w_idx < G->begin[v+1]; w_idx ++) {
	  node_t w = G->node_idx [w_idx];
	  if (isNeighbour(G,w,u)) {
	    data.coeff[v] += 1;
	  }
	  if (isNeighbour(G,u,w)) {
	    data.coeff[v] += 1;
	  }
	}
      }
      int neighbours = (int) (G->begin[v+1] - G->begin[v]);
      if(neighbours > 1) {
	data.coeff[v] = data.coeff[v]/(neighbours * (neighbours -1));
      }
    }
  }

  double clusterCoeff = 0;
  node_t v;
  for(v =0;v<G->numNodes;v++) {
    clusterCoeff += localClustering[v];
  }
  clusterCoeff = clusterCoeff/G->numNodes;
  data.clustercoeff = clusterCoeff;
  
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


graph* erdosRenyiGenerator(const GraphProperty& p) {
  int *stream1, *stream2;
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  edge_t maxEdges = 0;
  if(p.selfloop == false) {
    /* Should cover all edges 90 percent of time*/
    maxEdges = (edge_t) ( 1.1 * ((p.edgeProbablity) * p.numNodes) * (p.numNodes -1));
  } else {
    maxEdges = (edge_t) ( 1.1 * ((p.edgeProbablity) * p.numNodes) * (p.numNodes));
  }
  graph *G = allocatememoryforGraph(n, maxEdges);
   
  /* Write the no. of edges later */

  edge_t edg = 0;
  for (node_t i=0; i<p.numNodes; i++) {
    G->begin[i] = edg;
    for (node_t j=0; j<p.numNodes; j++) {
      if ((i==j) && (p.selfloop == false))		
	continue;
      if (p.edgeProbablity > sprng(stream1)) {
	if(weighted) {
	  int w = (int)( p.minWeight + ((double )(p.maxWeight - p.minWeight)) * sprng(stream2));
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[numNodes] = edg;
  assert(edg <= maxEdges);
  return G;
}

void choosePartition(node_t* u, node_t* v, node_t step, int* stream1) {
  double p = sprng(stream1);
  if (p < a) {
    
    /* Do nothing */
    
  } else if ((a < p) && (p < a+b)) {
    
    *v = *v + step;
    
  } else if ((a+b < p) && (p < a+b+c)) {
    
    *u = *u + step;
    
  } else if ((a+b+c < p) && (p < a+b+c+d)) {
    
    *u = *u + step;
    *v = *v + step;
  }

}
void varyParams(double* a, double* b, double* c, double* d, int* stream3,
		int* stream4) {
  double v, S;
  /* Allow a max. of 5% variation */
  v = 0.05;
  if (sprng(stream4) > 0.5)
    *a += *a * v * sprng(stream3);
  else 
    *a -= *a * v * sprng(stream3);
  
  if (sprng(stream4) > 0.5)
    *b += *b * v * sprng(stream3);
  else 
    *b += *b * v * sprng(stream3);
  
  if (sprng(stream4) > 0.5)
    *c += *c * v * sprng(stream3);
  else 
    *c -= *c * v * sprng(stream3);
  
  if (sprng(stream4) > 0.5)
    *d += *d * v * sprng(stream3);
  else 
    *d -= *d * v * sprng(stream3);
  // Normalize
  S = *a + *b + *c + *d;
  *a = *a/S;
  *b = *b/S;
  *c = *c/S;
  *d = *d/S;
  

}



graph* rmatGenerator(const GraphProperty& p) {
  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  int *stream1, *stream2, *stream3, *stream4;
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);
  stream3 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED3, SPRNG_DEFAULT);
  stream4 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED4, SPRNG_DEFAULT);
  double a0, b0, c0, d0;
  node_t u, v;
  for(edge_t it = 0; it< numEdges; it++) {
    a0 = p.a; b0 = p.b; c0 = p.c; d0= p.d;
    u = 0;
    v = 0;
    step = numNodes/2;
    while (step >= 1) {
      choosePartition(&u, &v, step, stream1);
      step = step / 2;
      varyParams(&a0, &b0, &c0, &d0, stream3, stream4);
    }

    if(p.selfloop == false && (u == v)) {
      it --;
      continue;
    }

    if(weighted) {
      w = (int)( p.minWeight + ((double )(p.maxWeight - minWeight)) * sprng(stream2));
      G->weights[edg] = w;
    }
    G->node_idx[i] = v;
    G->r_node_idx[i] = u;
    G->begin[i+1]++; // array start with 0
    
  }
  // sort and Normalize
  for(node_t n = 1;n< p.numNodes; n++) {
    G->begin[n+1] += G->begin[n];
  }

  doubleMergeSort(G->r_node_idx, G->node_idx, 0, p.numEdges);
  return G;  

}

graph* propertyControlledGraphGenerator(GraphProperty& p) {
  p.numNodes = (nodes_t) ( p.density / p.degree );
  if(densityflag == true && degreeflag == true) {
    p.numNodes = (nodes_t) ( p.density / p.degree );
    p.numEdges = (edges_t) ( p.degree * p.numNodes );
  } else if (densityflag == true) {
    p.numEdges = (edges_t) ( p.density * p.numNodes * p.numNodes);
  } else if (degreeflag == true) {
    p.numEdges = (edges_t) ( p.degree * p.numNodes);
  }
  // Associate expected clustering coefficient.
  // degree sd
  // aed.
  

  // TODO
  
}


/***
 * We can use this to scale up the community probabilities  
 * in the probability vector. Will be helpful when we are trying 
 * to generate graphs with large nmber of communities
 ***/
#define PVECTORSCALE 1

node_t getCommunityId(double* vec, double randComm, node_t len) {
  // Binary search tree
  node_t id;
  node_t base = 1;
  do {
    id = (len + base)/2 ;
    if(randComm < vec[id] && randComm > vec[id-1])
      return id;
    else if (randComm < vec[id])
      len = id;
    else
      base = id;
  } while(len > base);

  /** Debug statement can remove after enough testing **/
  assert(0 == 1);
  return NIL_NODE;
}



/***
 * probabilityVector: vector [length k+1, where k is the number of communities] is probability vector 
 * of a node to be part of community. k is the number fo communities  
 * edgeProbablilityMatrix: [size k*k] where each entry (i,j) gves the probability of edge between nodes of two
 * communities i and j.
 * k : number of communities
 **/
graph* SBMGraphGenerator(double* probablilityVector, double *edgeProbablilityMatrix, node_t k, const GraphProperty& p) {
  node_t* comm  = (node_t*) malloc (p.numNodes * sizeof(nodes_t));

  int* stream1,* stream2;
  
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  node_t i;
  for(i=0; i< numNodes; i++) {
    double randComm = sprng(stream1) * PVECTORSCALE;
    node_t commID = getCommunityId(probabilityVector, randComm, k);
    comm[i] = commID;
  }

  node_t j;
  /****
   * The expected number of edges is assumed to be saved in numEdges.
   * before this point.
   ****/
  maxEdges = (edges_t) (1.1 *p.numEdges);
  graph* G =  allocateMemoryforGraph(p.numNodes, maxEdges);
  node_t edg = 0;
  for(i=0; i< p.numNodes; i++) {
    G->begin[i] = edg;
    for (j=0; j< p.numNodes; j++) {
      if(p.selfloop == false && i == j)
	continue;
      double p = edgeProbabilityMatrix[(comm[i])*k + comm[j]];
      if (p > sprng(stream2)) {
	if(weighted) {
	  int w = (int)( minWeight + ((double)(maxWeight - minWeight)) * sprng(stream2));
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[p.numNodes] = edg;
  assert(maxEdges >= edg);
  free(probablilityVector);
  free(edgeProbablilityMatrix);
  return G;
}

graph* StocasticBlockModel() {
  
  // TODO
  
  
}

/****
 * k number of clusters
 * A is the intracluster edge probability
 * B is the intercluster edge probability
 * Ideally, A > B 
 */
graph* symmetricStocasticBlockModel(node_t k, double A, double B, GraphPorperty& p) {
  double *probvec = (double *) malloc((k+1) * sizeof(double));
  double *edgeProbMatrix = (double *) malloc ((k*k) * sizeof(double));
  double nodeprob = (double ) ( (double)PVECTORSCALE / k);


  double numClusters = (double) p.numNodes / k;
  /*** Calculate the expected number of edges ***/
  if(p.selfloop == true) {
    edge_t intraClusterEdges = (edge_t) numClusters * (k*k) * A;
    edge_t interClusterEdges = [(double) (p.numNodes * p.numNodes) - p.numClusters * (k*k)] * B;
    numEdges = intraClusterEdges + interClusterEdges;  
  } else {
    edge_t intraClusterEdges = (edge_t) numClusters * (k*(k-1)) * A;
    edge_t interClusterEdges = [(double) (p.numNodes * (p.numNodes - 1)) - p.numClusters * (k*(k -1))] * B;
    numEdges = intraClusterEdges + interClusterEdges;  
  }
  
  
  

#pragma omp parallel for
  for(node_t i=0;i<k;i++) {
    probvec[i] = nodeprob *i;
    node_t j;
    node_t base = i *k;
    for(j = 0; j<k; j++) {
      if(i == j)
	edgeProbMatrix[base + j] = A;
      else
	edgeProbMatrix[base +j ] = B
    }
  }
  probvec[k] = PVECTORSCALE;
  return SBMGraphGenerator(probvec, edgeProbMatrix, k);
}


graph* lowdegreeNetworkgraph(const GraphProperty& p) {

  int flag = 0;

  double scale = 0;
  nodes_t nNodes = (nodes_t) (p.degree/ p.density);
  edges_t edges = (edges_t) ((double)nNodes * p.degree);
  double vVal = log(nNodes/2);
  // (pi^2/6)
  if(p.degree < 3.28) {
    flag = 2;
    scale = p.degree/3.28;
  } else if(p.degree < vVal) {
    flag = 1;
    scale = p.degree/vVal;
  }

  // TODO
  assert(flag !=0);
  
  edge_t maxEdges = 0;
  
  assert(p.selfloop == 0);
  maxEdges = (edge_t) ( 1.1 * (edges));

  graph *G = allocatememoryforGraph(n, maxEdges);

  int *stream1, *stream2;
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  edge_t edg = 0;
  node_t aed;
  for (node_t i=0; i<p.numNodes; i++) {
    G->begin[i] = edg;
    for (node_t j=0; j<p.numNodes; j++) {
      if ((i==j))		
	continue;
      double prob = 0;
      if(i<j)
	aed = (i-j);
      else
	aed = (j-i);
      if(flag == 1) {
	prob =  scale * 1/aed;
      } else if(flag == 2){
	prob = scale / (aed*aed);
      }
      
      if (prob > sprng(stream1)) {
	if(weighted) {
	  int w = (int)( p.minWeight + ((double )(p.maxWeight - p.minWeight)) * sprng(stream2));
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[numNodes] = edg;
  assert(edg <= maxEdges);
  return G;    
}

graph* diagonalGraphGenerator(boolean randomizedmidpoint,const  GraphProperty& p) {
  edges_t numEdges = (edges_t) p.degree * p.numNodes;
  node_t i;
  maxEdges = (edges_t) (1.1 *numEdges);
  graph* G= allocateMemoryforGraph(p.numNodes, maxNodes);
  data.z1 = 0;
  data.generate = false;
  int* stream;
  stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);


  // Thread local data. 
  tldBoxMuller data;

  
  node_t edg = 0;
  for(i=0; i< p.numNodes; i++) {
    G->begin[i] = edg;
    node_t idegree = (node_t) generateGaussiandistriution(p.degree, p.degreesd, data);
    node_t midPoint; node_t startPoint;
    if(randomizedmidpoint) {
      midPoint = (node_t) ispring(stream1) % idegree;
    } else {
      midPoint = idegree/2;
    }

    if(midpoint < i) {
      midPoint = i;
    }
    if((idegree -midPoint) > (numNodes -i)) {
      midPoint = idgree - (numNodes -i);
    }
    startPoint = i - midPoint;
    node_t endPoint = startPoint + idegree;
    if(p.selfloop ==  false)
      endPoint++;
    for(node_t j =startPoint; j < endPoint;j++) {
      if ((i==j) && (p.selfloop == false))		
	continue;
      if(weighted) {
	int w = (int)( minWeight + ((double )(maxWeight - minWeight)) * sprng(stream2));
	G->weights[edg] = w;
      }
      G->node_idx[edg] = j;
      edg ++;
    }
  }
  return G;
}

//
// TODO have a global minima flag.
// if set we will accept random negative
// changes to improve the chances of avoiding
// falling into a local minima.
boolean globalMinimaFlag = false;
#define biasedCoinFlip 30
 
/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* adjustAED(graph *G, gdata& data, int* stream, int itrs, int flag) {
  // TODO
  data.changes = 0;
  data.negativeChanges = 0;
  int i;
  stream = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  j=0;
  node_t src[2*itrs];
  for(i=0; i< itrs; i++) {
    src[i+j] = (node_t) isprng(stream1) % numNodes;
    j++;
    src[i+j] = (node_t) isprng(stream1) % numNodes;
  }
  /***
   ** To parallelize these loops value dependencies shpuld be handled.
   **/
  node_t src1, src2;
  node_t buf[data.maxDegree];
  int bufWeight[data.maxDegree];
  for(int i =0; i< itrs; i++) {
    // swap i and itrs+i elements
    // if useful
    src1 =  src[i];
    src2 = src[itrs+i];
    node_t ed1, ed2;
    // calculate base ed
    node_t newed1 = 0;
    for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
      node_t d = G->node_idx[y];
      ed1 += src1 - d;
    }

    for(edge_t y= G->begin[src2];y < G->begin[src2+1]; y++) {
      node_t d = G->node_idx[y];
      ed2 += src2 - d;
    }
    
    // calculate potential new ed
    node_t newed1 = 0;
    for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
      node_t d = G->node_idx[y];
      newed1 += src2 - d;
    }

    node_t newed2 = 0;
    for(edge_t y= G->begin[src2];y < G->begin[src2+1]; y++) {
      node_t d = G->node_idx[y];
      newed2 += src1 - d;
    }
    bool swap = false;
    if( (((newed1 + newed2) >  (ed1+ ed2)) && (flag == 1)) ||
	(((newed1 + newed2) <  (ed1+ ed2)) && (flag == 0))) {
      swap  = true;
    } else if(globalMinimaFLag == true) {
      // swap randomly for negative change to avoid falling into
      // local minima
      if( isprng(stream1) % 100 < biasedCoinFlip) {
	swap = true;
	data.negativeChanges++
      }
    }

  
    if(swap == true && (src1 != src2)) {
      data.changes++;
      if( (G->begin[src1+1]-G->begin[src1]) < (G->begin[src2+1]-G->begin[src2])) {
	// shift right
	edge_t id =0;

	// copy src 2 edges to buf
	edge_t bufsize =  G->begin[src2+1] - G->begin[src2];
	for (edge_t y= G->begin[src2];y < G->begin[src2+1]; y ++) {
	  buf[id] = G->node_idx[y];
	  if(weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
	edge_t diff = (G->begin[src1+1]-G->begin[src1]) - (G->begin[src2+1]-G->begin[src2]);

	
	// move src1 edges to src2
	id = G->begin[src2  + diff];
	for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
	  G->node_idx[id] = G->node_idx[y];
	  if(weighted)
	    G->weights[id] = G->weights[y]; 
	  id++;
	}


	assert(id == G->begin[src2+1]);

	// shift the edges of nodes in between
	edge_t y = G->begin[src2 -1];
	for(; y > G->begin[src+1]; y--) {
	  G->node_idx[y+diff] = G->node_idx[y];
	  if(weighted)
	    G->weights[y+diff] = G->weights[y];
	}

	
	// finally copy back src2 edges to the place of src1
	id = G->begin[src1];
	for(edge_t y = 0; y< bufSize; y++) {
	  G->begin[id] = buf[y];
	  if(weighted)
	    G->weights[id] = bufWeight[y];
	  id++;
	}

	
	assert(id == (G->begin[src1+1] + diff));
	// Update begin values
	for(node_t n = src+1; n <= src2; n++) {
	  G->begin[n] += diff;
	}
    
	
	
      } else if((G->begin[src1+1]-G->begin[src1]) > (G->begin[src2+1]-G->begin[src2])) {
	// shift left
	edge_t id =0;
	// copy src 1 edges to buf
	edge_t bufsize =  G->begin[src1+1] - G->begin[src1];
	for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
	  buf[id] = G->node_idx[y];
	  if(weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
       	
	// move src2 edges to src1
	id = G->begin[src1];
	for (edge_t y= G->begin[src2];y < G->begin[src2+1]; y ++) {
	  G->node_idx[id] = G->node_idx[y];
	  if(weighted)
	    G->weights[id] = G->weights[y];
	  id++;
	}

	edge_t diff = (G->begin[src1+1]-G->begin[src1]) - (G->begin[src2+1]-G->begin[src2]);

	assert(id == G->begin[src1+1] - diff);

	// shift the edges of nodes in between
	for(node_t n = src1+1; n < src2; n++) {
	  for(edge_t y = G->begin[n]; y < G->begin[n+1];y++) {
	    G->node_idx[y-diff] = G->node_idx[y];
	    if(weighted)
	      G->weights[y-diff] = G->weights[y];
	  }
	}


	// finally copy back src1 edges to the place of src2
	id = G->begin[src2] - diff;
	for(edge_t y = 0; y< bufSize; y++) {
	  G->begin[id] = buf[y];
	  if(weighted)
	    G->weights[id] = bufWeight[y];
	  id++;
	}


	// Update begin values
	for(node_t n = src+1; n <= src2; n++) {
	  G->begin[n] -= diff;
	}
	
      } else {
	// equal degree
	// just swap
	edge_t id =0;
	
	// copy src 1 edges to buf
	edge_t bufsize =  G->begin[src1+1] - G->begin[src1];
	for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
	  buf[id] = G->node_idx[y];
	  if(weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
	
	
	// move src2 edges to src1
	id = G->begin[src1];
	for (edge_t y= G->begin[src2];y < G->begin[src2+1]; y ++) {
	  G->begin[id] = G->node_idx[y];
	  if(weighted)
	    G->weights[id] = G->weights[y];
	  id++;
	}
	// 
       assert(G->begin[src1+1]-G->begin[src1]) == (G->begin[src2+1]-G->begin[src2]);

	// finally copy back src1 edges to the place of src2
	id = G->begin[src2];
	for(edge_t y = 0; y< bufSize; y++) {
	  G->begin[id] = buf[y];
	  if(weighted)
	    G->weights[id] = bufWeight[y];
	  id++;
	}

	// No need to update begin values as
	// the degrees are same.
	/* // Update begin values */
	/* for(node_t n = src+1; n <= src2; n++) { */
	/*   G->begin[n] -= diff; */
	/* } */
      }


      
      //src1 < src2
      // assert(G->node_idx_src == NULL);
      // assert(G->r_begin ==  NULL);
      node_t id = 0;
      
      
      
    }

    data.aed + = ((double) (newed1 + newed2) -  (ed1 + ed2) )/ (numNodes * numEdges);
  }
  
  return G;
}


/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* changeClustering(graph* G, gdata& data, int* stream, int itrs, int flag) {
  // TODO increase intercluster edges.
}



/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* changeDegreesd(graph* G, gdata& data, int* stream, int itrs, int flag) {

  data.changes = 0;
  data.negativeChanges = 0;
  int i;
  stream = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  j=0;
  node_t src[2*itrs];
  for(i=0; i< itrs; i++) {
    src[i+j] = (node_t) isprng(stream1) % numNodes;
    j++;
    src[i+j] = (node_t) isprng(stream1) % numNodes;
  }

  // TODO
}
