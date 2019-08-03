#ifndef GRAPH_PROPERTY_H
#define GRAPH_PROPERTY_H

#include<string>

// TODO generate these values using scripts
typedef enum GraphModel {
  Random,
  ErdosRenyi,
  RMAT,
  SBM, // Stochastic Block Model
  SSBM, // Symmetric Stochastic Block Model
  LSSBM, // Log Sum Stochastic Block Model
  PCGG, // Property Controlled Graph Generation
  SSCA // GT graph generator. http://www.cse.psu.edu/~kxm85/software/GTgraph/gen.pdf
} GraphModel;


struct RMATdata{
  /***** RMAT **********/
  double a,b,c,d;
};




typedef enum GraphPropertyFlag {
  allproperty,
  degreeanddensity,
  clustercoeff,
  degreeaed,
  degreeanddegresd,
}GraphPropertyFlag;



class GraphProperty{
private:
  GraphModel model;
  GraphPropertyFlag property;
  bool weighted;
  node_t numNodes;
  edge_t numEdges;
  /*** Random Parameters ***/
  bool selfloop;
  /***** Erdos Renyi *****/
  double edgeProbability;
  /******** SSCA ********/
  node_t maxCliqueSize;
  node_t totClusters;
  /********** Graph prop value *******/
  /**
     NOTE: These parameters has higher preference over numEdges and numNodes.
     if these values contradict with the numEdges and numNodes 
     we update the numNodes and numEdges based on these values.
  **/
  double degree; // Average degree of each Node
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

  /****
   ** RMAT data
   ***/
  RMATdata rmat;
  
  bool degreeflag;
  bool densityflag;
  bool ccflag;
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

public:
  GraphProperty();
  bool updateConfigs(const std::string configfile);
  graph* callappropriategenerator();
  node_t get_numNodes(){return numNodes;}
  edge_t get_numEdges(){return numEdges;}
  int get_minWeight() {return minWeight;}
  int get_maxWeight(){return maxWeight;}
  double get_edgeProbability() {return edgeProbability;}
  double get_degree(){return degree;}
  double get_degreesd(){return degreesd;}
  double get_clusteringCoeff() {return clusteringCoeff;}
  bool get_selfloop(){return selfloop;}
  bool get_weighted(){return weighted;}
  double get_aed() {return aed;}
  double get_density() {return density;}
  RMATdata& get_RMATdata() {return rmat; }
  void set_numNodes(node_t n) {numNodes = n;}
};



// TO CHECK: double or float of individual degree, coeff and ed ?
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

#endif
