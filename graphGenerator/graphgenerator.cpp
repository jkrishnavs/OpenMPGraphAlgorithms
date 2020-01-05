/*****
 * Original base algorithm from GTGraph. 
 ****/
#include <unistd.h>
#include <string.h>
#include <sprng_cpp.h>

#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "graphProperty.hh"
#include "normaldistribution.h"
#include "randomgraph.h"
#include "graphgenerator.h"
#include "graphprop.h"


/**
 * Dummy function since the generator is linked to main functions.
 */
void kernel(graph *G) { /* empty*/}

class graph_model {
protected:
  GraphProperty * p; 
  GraphModel model;
public:
  virtual graph* generate_graph()= 0;
  virtual bool check_feasibility() = 0;
  graph_model(GraphProperty *prop) {
    p = prop;
  }
};

class RandomGenerator: public graph_model {
public:
  graph* generate_graph();
  RandomGenerator(GraphProperty *p): graph_model(p) {
    model = GraphModel::Random;
  }
  bool check_feasibility(){return true;}
};


class ErdosRenyiGenerator: public graph_model {
public:
  graph* generate_graph();
  ErdosRenyiGenerator(GraphProperty *p): graph_model(p) {
    model = GraphModel::ErdosRenyi;
  }
  bool check_feasibility(){return true;}
};



class RMATGenerator: public graph_model {
public:
  graph* generate_graph();
  RMATGenerator(GraphProperty *p): graph_model(p){
    model = GraphModel::RMAT;
  }
  bool check_feasibility(){return true;}

};


class SBMGenerator: public graph_model {
public:
  graph* generate_graph();
  SBMGenerator(GraphProperty* p):graph_model(p) {
    model = GraphModel::SBM;
  }
  bool check_feasibility(){return true;}
  graph* SBMGraphGenerator(double* probabilityVector, double *edgeProbabilityMatrix, node_t k);
};

class SSBMGenerator: public SBMGenerator {
public:
  graph* generate_graph();
  SSBMGenerator(GraphProperty* p):SBMGenerator(p) {
    model = GraphModel::SSBM;
  }
  edge_t expected_edges(bool selfloop, node_t numClusters, double k, double A);
  bool check_feasibility();
};

class LowDegreeNetworkGraph: public graph_model{
 public:
  LowDegreeNetworkGraph(GraphProperty *p): graph_model(p){
    model =  GraphModel::LDNG;
  }
  graph *generate_graph(); 
  bool check_feasibility() {return true;}
};

class DiagonalGraphGenerator: public graph_model{
private:
  bool randomizedmidpoint;
public:
  DiagonalGraphGenerator(GraphProperty *p): graph_model(p) {
    model = GraphModel::DGG;
    randomizedmidpoint = false;
  }
  graph *generate_graph();
  bool check_feasibility() {return true;}  
};

graph* callappropriategenerator(GraphProperty* p) {
  graph_model* gm = NULL;
  if(p->get_model() == GraphModel::Random) {
    gm = new RandomGenerator(p);
  } else if (p->get_model()  ==   GraphModel::ErdosRenyi) {
    gm = new ErdosRenyiGenerator(p);
  } else if (p->get_model() == GraphModel::RMAT) {
    gm = new RMATGenerator(p);
  } else if (p->get_model() == GraphModel::SSCA) {
    gm = NULL; //TODO
  } else if (p->get_model() == GraphModel::SBM ||
	     p->get_model() == GraphModel::LSSBM ||
	     p->get_model() == GraphModel::PCGG) {
    gm = new SBMGenerator(p);
  } else if (p->get_model() == GraphModel::SSBM) {
    gm = new SSBMGenerator(p);
  } else if (p->get_model() == GraphModel::DGG) {
    gm = new DiagonalGraphGenerator(p);
  } else if (p->get_model() == GraphModel::LDNG) {
    gm = new LowDegreeNetworkGraph(p);
  }
  bool feasible = gm->check_feasibility();
  if(feasible) {
    return gm->generate_graph();
  }
  return NULL;
}



  

/**
 * The main function
 */ 
int runalgo(int argc, char** argv) {
  if(argc != 4) {
    const char* argList[3] = {" <configfile.conf>" , "<outputfile.edges>","<propfile.prop>"};
    printError(INCORRECT_ARG_LIST, 3, argList);
    return -1;
  }
  GraphProperty* p = new GraphProperty();
  std::string configfile = argv[1]; 
  bool scanFlag = p->updateConfigs(configfile);
  if(scanFlag == false) {
    printError(INCORRECT_CONFIG_FILE, 3, NULL);
    return -1;
  }
  graph * G = callappropriategenerator(p);
  assert(G != NULL);
  
  writeBackGraph(G, argv[2]);
  double avgincident =  (double)G->numEdges / G->numNodes;
  /*** Graph prop collection ***/
  printf("Avg Adjecency =  %f \n", avgincident);
  printf("\nAvg ClusterCoeff = %f \n", clusterCoeff);
  printf("Avg Edge Distance = %f \n", aed);
  printf("Sparsity = %.7f \n", sparsityMeasure);
  printf("The Scaled percentange Triangles = %f\n\n", scaledT);

  return 0;
}


inline int generateweight(GraphProperty& p, Sprng* stream2) {
  return (int)( p.get_minWeight() + ((double )(p.get_maxWeight() - p.get_minWeight())) * stream2->sprng());
}



/*
graph* SSCAGenerator(GraphProperty& p) {

  return NULL;
}
*/


#define cmrg 3
graph* RandomGenerator::generate_graph() {
  FILE* outfp;
  
  Sprng* stream1 =  SelectType(cmrg);
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  Sprng* stream2 = SelectType(cmrg);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);
  graph* G = allocateMemoryforGraph(p->get_numNodes(), p->get_numEdges(), p->get_weighted());
  G->r_node_idx = (node_t*) malloc (p->get_numEdges() * sizeof(node_t));
  G->numNodes = p->get_numNodes();
  G->numEdges = p->get_numEdges();
  //  #pragma omp parallel for
  for(edge_t i = 0; i< p->get_numEdges(); i++) {
    node_t u = (node_t) stream1->isprng() % p->get_numNodes();
    node_t v = (node_t) stream1->isprng() % p->get_numNodes();
    if((u == v) && (p->get_selfloop() == false)) {
      i--;
      continue;
    }
    if(p->get_weighted()) {
      int w = generateweight(*p,stream2);
      G->weights[i] = w;
    }
    G->node_idx[i] = v;
    G->r_node_idx[i] = u;
    //#pragma omp atomic
    G->begin[u]++; // array start with 0
  }
  for(node_t n = 1;n< p->get_numNodes(); n++) {
    G->begin[n+1] += G->begin[n];
  }
  doubleMergeSort(G->r_node_idx, G->node_idx, 0, p->get_numEdges());
  free(G->r_node_idx);
  delete stream1;
  delete stream2;  
  return G;
}

#define compare_double_list(l1,l2,i1,i2) ((l1[i2]< l1[i1]) ? true: ( (l1[i2]== l1[i1]) ? l2[i2] < l2[i1]: false))

#define compare_DD_list(t1,t2, l1, l2, ti, i) ((l1[i]< t1[ti]) ? true: ( (l1[i]== t1[ti]) ? t2[ti] < l2[i]: false))

inline void swap_double_list(node_t* l1, node_t* l2, edge_t i1, edge_t i2) {
  node_t tmp;
  tmp = l1[i1]; l1[i1] = l1[i2]; l1[i2] = tmp;
  tmp = l2[i1]; l2[i1] = l2[i2]; l2[i2] = tmp;
}

inline void copy_double_list(node_t* s1, node_t* s2, node_t* d1, node_t* d2, edge_t si, edge_t di) {
  d1[di] = s1[si];
  d2[di] = s2[si];
}


/*** Update Graph Meta data **/
void updateGdata(graph* G, gdata& data) {
  data.changes = 0;
  data.negativeChanges = 0;
  data.maxDegree = 0;

  int maxVal = 0;
  long long degreeSum = 0;

#pragma omp parallel for reduction(+: degreeSum) reduction(max: maxVal)
  for (node_t v = 0; v < G->numNodes; v ++) {
    if((G->begin[v+1] - G->begin[v]) > maxVal) {
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
  for (node_t v = 0; v < G->numNodes; v ++) {
    xval = (G->begin[v+1] - G->begin[v]);
    degreesd +=  (xval - data.degree) * (xval -data.degree) ;
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

  double *localClustering = (double*)malloc (G->numNodes * sizeof(double));
  // TODO calculate Local Clustering
  double clusterCoeff = 0;
  node_t v;
  for(v =0;v<G->numNodes;v++) {
    clusterCoeff += localClustering[v];
  }

  free(localClustering);
  clusterCoeff = clusterCoeff/G->numNodes;
  data.clustercoeff = clusterCoeff;
  
}

#define _SELECT_LOWER 0
#define _SELECT_UPPER 1
#define _SELECT_TEMP 2 


void merge(node_t* l1, node_t* l2, edge_t start, edge_t mid, edge_t end) {
  node_t t1 = start;  node_t t2 = mid;
  node_t tmp;  node_t temp; node_t w = mid-1;
  /*
   * We are allocating some buffer to copy out the first half of the 
   * array (start- mid).
   * this buffer acts as a queue, FIFO.
   * with its front pointed to by tpf and back pointed to by 
   * tpb. i.e. the queue is empty when tpf == tpb. 
   */
  node_t tempqueue[(mid - start)*2]; // should we use malloc ?
  node_t tpf = 0; node_t tpb = 0; 
  node_t offset = mid-start;
  node_t *tq2 = &(tempqueue[offset]);

  
  int _flag;
  while(t2< end) {
    if(tpf == tpb && t1 < mid) {
      if(compare_double_list(l1, l2, t1, t2)) {
	_flag = _SELECT_UPPER;
      } else {
	_flag = _SELECT_LOWER;
      }
    } else{
      if(compare_DD_list(tempqueue, tq2, l1, l2, tpb, t2)) {
	_flag = _SELECT_UPPER;
      } else {
	_flag = _SELECT_TEMP;
      }
    }
    // If upper or temp copy the list to current positions
    if(t1 < mid && (_flag == _SELECT_UPPER || _flag == _SELECT_TEMP)) {
      // backup lower to temp
      	copy_double_list(l1, l2, tempqueue, tq2, t1, tpf);
	tpf++;
    }
    if(_flag == _SELECT_UPPER) {
	copy_double_list(l1, l2, l1, l2, t2, t1);
	t2++;
    } else if (_flag == _SELECT_TEMP) {
	copy_double_list(tempqueue, tq2, l1, l2,  tpb, t1);
	tpb++;
    }
    if(t1> mid  && tpf == tpb) {
      // the loop is already sorted.
      t1 = t2; // for the assert.
      break;
    }
    t1++;
  }
  if(t1 < mid) {
    assert(0 ==1);
  }
  // Just copy back the buffered lower half to the end of the merged list.
  while(tpf > tpb ) {
    copy_double_list(tempqueue, tq2, l1, l2,  tpb, t1);
    tpb++;
    t1++;
  }
 
  assert(tpf == tpb);
  assert(t1 == t2); 
}

void doubleMergeSort(node_t* l1, node_t* l2, edge_t left, edge_t right) {
  if((right - left)> 10) {
    edge_t mid = left + (right-left)/2;
    doubleMergeSort(l1,l2, left, mid);
    doubleMergeSort(l1,l2, mid, right);
    merge( l1, l2, left, mid, right);
  } else {
    edge_t itr = left;
    for(; itr < right; itr++) {
      edge_t itr2= itr;
      for(; itr2 < right; itr2++) {
	if(compare_double_list(l1,l2, itr, itr2)) {
	  swap_double_list(l1,l2, itr, itr2);
	}
      }
    }
  }
}

graph* allocateMemoryforGraph(node_t n, edge_t m, bool weighted) {
  graph *G = (graph*) malloc (sizeof(graph));
  G->begin = (edge_t*) malloc( (n+1)* sizeof(edge_t));
  G->node_idx = (node_t*) malloc(m * sizeof(node_t));
  if(weighted == true) {
    G->weighted = true;
    G->weights = (int*) malloc(m* sizeof(int));
  }
  /* Initialize */
  memset(G->begin, 0, (n+1) * sizeof(edge_t));
  return G;
}


graph* ErdosRenyiGenerator::generate_graph() {
  Sprng *stream1 = SelectType(cmrg);
  Sprng *stream2 = SelectType(cmrg);
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);
  if(p->get_edgeProbability() != 0) {
    if(p->get_selfloop() == false){
      p->set_numEdges( (edge_t)(p->get_edgeProbability() * p->get_numNodes()) * (p->get_numNodes() -1));
    } else {
      p->set_numEdges( (edge_t)(p->get_edgeProbability() * p->get_numNodes()) * p->get_numNodes());
    }
  } else {
    assert(p->get_numEdges() != 0);
    if(p->get_selfloop() == false){
      p->set_edgeProbability((double)(p->get_numEdges())/ (p->get_numNodes() * (p->get_numNodes() -1)));
    } else {
      p->set_edgeProbability((double)(p->get_numEdges())/ (p->get_numNodes() * p->get_numNodes()));
    }
  }

  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  edge_t maxEdges = 1.1 * p->get_numEdges();
  graph *G = allocateMemoryforGraph(p->get_numNodes(),
				    maxEdges, p->get_weighted());
   
  /* Write the no. of edges later */
  node_t numNodes = p->get_numNodes();
  G->numNodes  = numNodes;
  edge_t edg = 0;
  for (node_t i=0; i<numNodes; i++) {
    G->begin[i] = edg;
    for (node_t j=0; j<numNodes; j++) {
      if ((i==j) && (p->get_selfloop() == false))		
	continue;
      if (p->get_edgeProbability() > stream1->sprng()) {
	if(p->get_weighted()) {
	  int w = generateweight(*p,stream2);
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[numNodes] = edg;
  G->numEdges = edg;
  assert(edg <= maxEdges);
  delete stream1;
  delete stream2;
  return G;
}

void choosePartition(RMATdata& r, node_t* u, node_t* v, node_t step, Sprng* stream1) {
  double probability = stream1->sprng();
  if (probability <=r.a) {
    
    /* Do nothing */
    
  } else if (probability <= (r.a+r.b)) {
    
    *v = *v + step;
    
  } else if (probability <= (r.a+r.b + r.c)) {
    
    *u = *u + step;
    
  } else {
    
    *u = *u + step;
    *v = *v + step;
  }

}
void varyParams(double* a, double* b, double* c, double* d,
		Sprng* stream3, Sprng* stream4) {
  double v, S;
  /* Allow a max. of 5% variation */
  v = 0.05;
  if (stream4->sprng() > 0.5)
    *a += *a * v * stream3->sprng();
  else 
    *a -= *a * v * stream3->sprng();
  
  if (stream4->sprng() > 0.5)
    *b += *b * v * stream3->sprng();
  else 
    *b += *b * v * stream3->sprng();
  
  if (stream4->sprng() > 0.5)
    *c += *c * v * stream3->sprng();
  else 
    *c -= *c * v * stream3->sprng();
  
  if (stream4->sprng() > 0.5)
    *d += *d * v * stream3->sprng();
  else 
    *d -= *d * v * stream3->sprng();
  // Normalize
  S = *a + *b + *c + *d;
  *a = *a/S;
  *b = *b/S;
  *c = *c/S;
  *d = *d/S;
  

}



graph* RMATGenerator::generate_graph() {
  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  Sprng *stream1 = SelectType(cmrg);
  Sprng *stream2 = SelectType(cmrg);
  Sprng *stream3 = SelectType(cmrg);
  Sprng *stream4 = SelectType(cmrg);
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);
  stream3->init_sprng(0, 1, SPRNG_SEED3, SPRNG_DEFAULT);
  stream4->init_sprng(0, 1, SPRNG_SEED4, SPRNG_DEFAULT);
  double a0, b0, c0, d0;
  node_t u, v;
  RMATdata  rmat = p->get_RMATdata();
  graph* G = allocateMemoryforGraph(p->get_numNodes(), p->get_numEdges(), p->get_weighted());
  /* allocate memory to store the edge sources*/
  G->r_node_idx = (node_t*) malloc (p->get_numEdges()* sizeof(node_t));


  for(edge_t it = 0; it< p->get_numEdges(); it++) {
    a0 = rmat.a; b0 = rmat.b; c0 = rmat.c; d0= rmat.d;
    u = 0;
    v = 0;
    node_t step = p->get_numNodes()/2;
    while (step >= 1) {
      choosePartition(rmat, &u, &v, step, stream1);
      step = step / 2;
      varyParams(&a0, &b0, &c0, &d0, stream3, stream4);
    }

    if(p->get_selfloop() == false && (u == v)) {
      it --;
      continue;
    }

    if(p->get_weighted()) {
      int w = generateweight(*p, stream2);
      G->weights[it] = w;
    }
    G->node_idx[it] = v;
    G->r_node_idx[it] = u;
    G->begin[u]++;
  }
  // sort and normalize
  for(node_t n = 1;n< p->get_numNodes(); n++) {
    G->begin[n+1] += G->begin[n];
  }
  G->numNodes  = p->get_numNodes();
  G->numEdges = p->get_numEdges();
  doubleMergeSort(G->r_node_idx, G->node_idx, 0, p->get_numEdges());
  free(G->r_node_idx);
  delete stream1;
  delete stream2;
  delete stream3;
  delete stream4;
  return G;  

}

/* graph* propertyControlledGraphGenerator(GraphProperty& p) { */
/*   p->get_numNodes = (node_t) ( p->get_density() / p->get_degree() ); */
/*   if(densityflag == true && degreeflag == true) { */
/*     p->get_numNodes = (node_t) ( p->get_density{} / p->get_degree() ); */
/*     p->get_numEdges() = (edge_t) ( p->get_degree() * p->get_numNodes ); */
/*   } else if (densityflag == true) { */
/*     p->get_numEdges() = (edge_t) ( p->get_density() * p->get_numNodes * p->get_numNodes); */
/*   } else if (degreeflag == true) { */
/*     p->get_numEdges() = (edge_t) ( p->get_degree() * p->get_numNodes); */
/*   } */
/*   // Associate expected clustering coefficient. */
/*   // degree sd */
/*   // aed. */
  

  
/* } */


/***
 * We can use this to scale up the community probabilities  
 * in the probability vector. Will be helpful when we are trying 
 * to generate graphs with large number of communities
 ***/
#define PVECTORSCALE 1

node_t getCommunityId(double* vec, double randComm, node_t len) {
  node_t id;
  node_t base = 0;
  do {
    id = (len + base)/2 ;
    if(randComm < vec[id+1] && randComm > vec[id])
      return id;
    else if (randComm < vec[id])
      len = id;
    else
      base = id+1;
  } while(len > base);
  /** Debug statement can remove after enough testing **/
  assert(0 == 1);
  return NIL_NODE;
}



/***
 * probabilityVector: vector [length k+1, where k is the number of communities]
 * is probability vector 
 * of a node to be part of community. k is the number of  communities. 
 * To speed up the process  the values of 
 * the vector is left tail appended. 
 * edgeProbabilityMatrix: [size k*k] where each entry (i,j) gves the probability of edge between nodes of two
 * communities i and j.
 * k : number of communities
 **/
graph* SBMGenerator::SBMGraphGenerator(double* probabilityVector, double *edgeProbabilityMatrix, node_t k) {
  node_t* comm  = (node_t*) malloc (p->get_numNodes() * sizeof(node_t));

  

  Sprng* stream1 = SelectType(cmrg);
  Sprng* stream2 = SelectType(cmrg);
  
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  node_t i;

  for(i=0; i< p->get_numNodes(); i++) {
    double randComm = stream1->sprng() * PVECTORSCALE;
    node_t commID = getCommunityId(probabilityVector, randComm, k+1);
    comm[i] = commID;
  }
  

  node_t j;
  /****
   * The expected number of edges is assumed to be saved in numEdges.
   * before this point.
   ****/
  edge_t maxEdges = (edge_t) (1.1 *p->get_numEdges());
  graph* G =  allocateMemoryforGraph(p->get_numNodes(), maxEdges, p->get_weighted());
  node_t edg = 0;
  for(i=0; i< p->get_numNodes(); i++) {
    G->begin[i] = edg;
    for (j=0; j< p->get_numNodes(); j++) {
      if(p->get_selfloop() == false && i == j)
	continue;
      int idx = (comm[i])*k + comm[j];
      double probability = edgeProbabilityMatrix[idx];
      double generated_value = stream2->sprng();
      if (probability > generated_value) {
	if(p->get_weighted()) {
	  int w = generateweight(*p,stream2);
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[p->get_numNodes()] = edg;
  G->numNodes = p->get_numNodes();
  G->numEdges = edg;
  assert(maxEdges >= edg);
  free(probabilityVector);
  free(edgeProbabilityMatrix);
  //delete stream1;
  //delete stream2;
  return G;
}


graph* SBMGenerator::generate_graph() {


  return NULL;
}  
 //  if(p->get_degreesd() < 0.5) {
//     // symmetric model
//     // if clustering coefficent is high
//     // then the cluster size is closer to
//     //

//     node_t clusterSize = (node_t) p->get_degree()+1;
//     node_t numNodes = (node_t) ( p->get_degree()/ p->get_density());
//     node_t k = (node_t) (numNodes / clusterSize);
//     p->set_numNodes(numNodes);
//     double sizefrac = ((double) clusterSize / numNodes);

//     double edgprop = cbrt(p->get_clusteringCoeff() * 0.9);
//     // currently decided 10% of clustering coefficnt value should
//     // come from iter cluster edges.
//     // TODO make it flexible
//     // for more details on derivations see sbm.pdf in docs
//     // folder
//     double A = edgprop;
//     assert(0< A && A < 1);
//     double B = (1 -A) * (clusterSize/numNodes);
//     graph *G = symmetricStocasticBlockModel(k, A, B);
//     if(p->get_aed() < 0.2) {
//       // TODO sort
//     }
    
//   } else {
//     // TODO LSSBM,SBM, PCGG
//   }
  
// }

bool SSBMGenerator::check_feasibility() {
  edge_t intraClusterEdges, interClusterEdges, numEdges;
  SSBMdata ssbm = p->get_SSBMdata();
  node_t k = ssbm.k;
  double A = ssbm.alpha;
  double B = ssbm.beta;


  double sizeofClusters = (double) p->get_numNodes() / k;
  /*** Calculate the expected number of edges ***/
  
  if(p->get_selfloop() == true) {
    double clusterEdges = (sizeofClusters * sizeofClusters) * k; 
    intraClusterEdges =  clusterEdges * A;
    interClusterEdges = ((double) (p->get_numNodes() * p->get_numNodes()) - clusterEdges) * B;
    numEdges = intraClusterEdges + interClusterEdges;  
  } else {
    double clusterEdges = (sizeofClusters * (sizeofClusters-1)) * k; 
    intraClusterEdges =  clusterEdges * A;
    interClusterEdges = ((double) (p->get_numNodes() * (p->get_numNodes()-1)) - clusterEdges) * B;
    numEdges = intraClusterEdges + interClusterEdges;  
  }
  edge_t maxEdges = 1.05 * p->get_numEdges();
  //printf("Min %d  Max %d clsuters %lf k = %d \n", numEdges, maxEdges, sizeofClusters, k);
  if(numEdges  > maxEdges) {
    const char *msg[] = {"ALERT:We need to scale down the intra and inter cluster  ratio as the  number of edges that will be generated is much more than the  permissible limit."} ;
    printMsg(1, msg);
    double scale = (double) numEdges/ p->get_numEdges();
    p->update_SSBMdata(ssbm.alpha/scale, ssbm.beta/scale, ssbm.k);
  }
  return true;
}

/****
 * k number of clusters
 * A is the intracluster edge probability
 * B is the intercluster edge probability
 * Ideally, A > B 
 */
graph* SSBMGenerator::generate_graph() {
  SSBMdata ssbm = p->get_SSBMdata();
  node_t k = ssbm.k;
  double A = ssbm.alpha;
  double B = ssbm.beta;

  double *probvec = (double *) malloc((k+1) * sizeof(double));
  double *edgeProbMatrix = (double *) malloc ((k*k) * sizeof(double));
  double nodeprob =  ( (double)PVECTORSCALE / k);

  
#pragma omp parallel for
  for(node_t i=0;i<k;i++) {
    probvec[i] = nodeprob * i;
    node_t j;
    node_t base = i *k;
    for(j = 0; j<k; j++) {
      if(i == j)
	edgeProbMatrix[base + j] = A;
      else
	edgeProbMatrix[base +j ] = B;
    }
  }

  
  probvec[k] = PVECTORSCALE;
  return SBMGraphGenerator(probvec, edgeProbMatrix, k);
}


graph* LowDegreeNetworkGraph::generate_graph(){

  int flag = 0;
  double scale = 0;
  node_t numNodes = (node_t) (p->get_degree()/ p->get_density());
  edge_t edges = (edge_t) ((double)numNodes * p->get_degree());
  double vVal = log(numNodes/2);
  // (pi^2/6)
  if(p->get_degree() < 3.28) {
    flag = 2;
    scale = p->get_degree()/3.28;
  } else if(p->get_degree() < vVal) {
    flag = 1;
    scale = p->get_degree()/vVal;
  }
  
  assert(flag !=0);
  edge_t maxEdges = 0;
  assert(p->get_selfloop() == 0);
  maxEdges = (edge_t) ( 1.1 * (edges));

  graph *G = allocateMemoryforGraph(numNodes, maxEdges, p->get_weighted());

  Sprng *stream1 = SelectType(cmrg);
  Sprng *stream2 = SelectType(cmrg);
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

  edge_t edg = 0;
  node_t aed;
  for (node_t i=0; i<p->get_numNodes(); i++) {
    G->begin[i] = edg;
    for (node_t j=0; j<p->get_numNodes(); j++) {
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
      
      if (prob > stream1->sprng()) {
	if(p->get_weighted()) {
	  int w = generateweight(*p,stream2);
	  G->weights[edg] = w;
	}
	G->node_idx[edg] = j;
	edg ++;	  
      }
    }
  }
  G->begin[numNodes] = edg;
  assert(edg <= maxEdges);
  delete stream1;
  delete stream2;
  return G;    
}

graph*   DiagonalGraphGenerator::generate_graph(){
  edge_t numEdges = (edge_t) p->get_degree() * p->get_numNodes();
  node_t i;
  node_t numNodes = p->get_numNodes();
  edge_t maxEdges = (edge_t) (1.1 *p->get_numEdges());
  graph* G= allocateMemoryforGraph(p->get_numNodes(), maxEdges, p->get_weighted());
  normaldistribution::tldBoxMuller data;
  data.z1 = 0;
  data.generate = false;
  Sprng* stream1 = SelectType(cmrg);
  Sprng* stream2 = SelectType(cmrg);
  stream1->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  stream2->init_sprng(0, 1, SPRNG_SEED2, SPRNG_DEFAULT);


  
  
  node_t edg = 0;
  for(i=0; i< p->get_numNodes(); i++) {
    G->begin[i] = edg;
    node_t idegree = (node_t) normaldistribution::generateGaussiandistribution(p->get_degree(), p->get_degreesd(), data);
    node_t midPoint; node_t startPoint;
    if(randomizedmidpoint) {
      midPoint = (node_t) stream1->isprng() % idegree;
    } else {
      midPoint = idegree/2;
    }

    if(midPoint < i) {
      midPoint = i;
    }
    if((idegree -midPoint) > (numNodes -i)) {
      midPoint = idegree - (numNodes -i);
    }
    startPoint = i - midPoint;
    node_t endPoint = startPoint + idegree;
    if(p->get_selfloop() ==  false)
      endPoint++;
    for(node_t j =startPoint; j < endPoint;j++) {
      if ((i==j) && (p->get_selfloop() == false))		
	continue;
      if(p->get_weighted()) {
	int w = generateweight(*p, stream2);
	G->weights[edg] = w;
      }
      G->node_idx[edg] = j;
      edg ++;
    }
  }
  delete stream1;
  delete stream2;
  return G;
}

#define biasedCoinFlip 30
 
/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* adjustAED(graph *G, gdata& data, Sprng* stream, int itrs, int flag) {
//
// We have a global minima flag.
// if set we will accept random negative
// changes to improve the chances of avoiding
// falling into a local minima.
  bool globalMinimaFlag = false;
  data.changes = 0;
  data.negativeChanges = 0;
  int i;
  //stream->init_prng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  int j=0;
  node_t src[2*itrs];
  for(i=0; i< itrs; i++) {
    src[i+j] = (node_t) stream->isprng() % G->numNodes;
    j++;
    src[i+j] = (node_t) stream->isprng() % G->numNodes;
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
    } else if(globalMinimaFlag == true) {
      // swap randomly for negative change to avoid falling into
      // local minima
      if( stream->isprng() % 100 < biasedCoinFlip) {
	swap = true;
	data.negativeChanges++;
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
	  if(G->weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
	edge_t diff = (G->begin[src1+1]-G->begin[src1]) - (G->begin[src2+1]-G->begin[src2]);

	
	// move src1 edges to src2
	id = G->begin[src2  + diff];
	for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
	  G->node_idx[id] = G->node_idx[y];
	  if(G->weighted)
	    G->weights[id] = G->weights[y]; 
	  id++;
	}


	assert(id == G->begin[src2+1]);

	// shift the edges of nodes in between
	edge_t y = G->begin[src2 -1];
	for(; y > G->begin[src1+1]; y--) {
	  G->node_idx[y+diff] = G->node_idx[y];
	  if(G->weighted)
	    G->weights[y+diff] = G->weights[y];
	}

	
	// finally copy back src2 edges to the place of src1
	id = G->begin[src1];
	for(edge_t y = 0; y< bufsize; y++) {
	  G->begin[id] = buf[y];
	  if(G->weighted)
	    G->weights[id] = bufWeight[y];
	  id++;
	}

	
	assert(id == (G->begin[src1+1] + diff));
	// Update begin values
	for(node_t n = src1+1; n <= src2; n++) {
	  G->begin[n] += diff;
	}
    
	
	
      } else if(( (G->begin[src1+1]) - (G->begin[src1]) ) >
		( (G->begin[src2+1]) - (G->begin[src2]) )) {
	// shift left
	edge_t id =0;
	// copy src 1 edges to buf
	edge_t bufsize =  G->begin[src1+1] - G->begin[src1];
	for (edge_t y= G->begin[src1];y < G->begin[src1+1]; y ++) {
	  buf[id] = G->node_idx[y];
	  if(G->weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
       	
	// move src2 edges to src1
	id = G->begin[src1];
	for (edge_t y= G->begin[src2];y < G->begin[src2+1]; y ++) {
	  G->node_idx[id] = G->node_idx[y];
	  if(G->weighted)
	    G->weights[id] = G->weights[y];
	  id++;
	}

	edge_t diff = (G->begin[src1+1]-G->begin[src1]) - (G->begin[src2+1]-G->begin[src2]);

	assert(id == G->begin[src1+1] - diff);

	// shift the edges of nodes in between
	for(node_t n = src1+1; n < src2; n++) {
	  for(edge_t y = G->begin[n]; y < G->begin[n+1];y++) {
	    G->node_idx[y-diff] = G->node_idx[y];
	    if(G->weighted)
	      G->weights[y-diff] = G->weights[y];
	  }
	}


	// finally copy back src1 edges to the place of src2
	id = G->begin[src2] - diff;
	for(edge_t y = 0; y< bufsize; y++) {
	  G->begin[id] = buf[y];
	  if(G->weighted)
	    G->weights[id] = bufWeight[y];
	  id++;
	}


	// Update begin values
	for(node_t n = src1+1; n <= src2; n++) {
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
	  if(G->weighted)
	    bufWeight[id] = G->weights[y];
	  id++;
	}
	
	
	
	// move src2 edges to src1
	id = G->begin[src1];
	for (edge_t y= G->begin[src2];y < G->begin[src2+1]; y ++) {
	  G->begin[id] = G->node_idx[y];
	  if(G->weights != NULL)
	    G->weights[id] = G->weights[y];
	  id++;
	}
	// 
	assert( (G->begin[src1+1]-G->begin[src1]) == (G->begin[src2+1]-G->begin[src2]));

	// finally copy back src1 edges to the place of src2
	id = G->begin[src2];
	for(edge_t y = 0; y< bufsize; y++) {
	  G->begin[id] = buf[y];
	  if(G->weighted)
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

    data.aed += ((double) (newed1 + newed2) -  (ed1 + ed2) )/ (G->numNodes * G->numEdges);
  }
  
  return G;
}


/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* changeClustering(graph* G, gdata& data, Sprng* stream, int itrs, int flag) {
  
  // TODO increase intercluster edges.
}



/* Increase or decrese aed
 * if flag = 0 decrease
 * if flag = 1 increase
*/
graph* changeDegreesd(graph* G, gdata& data, Sprng* stream, int itrs, int flag) {

  /*
  data.changes = 0;
  data.negativeChanges = 0;
  int i;
  Sprng* stream->init_sprng(0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
  j=0;
  node_t src[2*itrs];
  for(i=0; i< itrs; i++) {
    src[i+j] = (node_t) stream->isprng() % G->numNodes;
    j++;
    src[i+j] = (node_t) stream->isprng() % G->numNodes;
  }

  delete
  */

  return G;
}
