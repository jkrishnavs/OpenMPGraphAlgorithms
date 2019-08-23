#include<unistd.h>
#include<cstring>
#include "graph.h"
#include "print.h"
#include "graphProperty.hh"
//#include "graphgenerator.h" 


/**** default configs ******/
GraphProperty::GraphProperty() {
  weighted = false; // wf
  selfloop = false; // sl
  sortedSBMGraph = true; 

  degreeflag = false;
  densityflag = false;
  ccflag = false;
  aedflag = false;
  degreesdflag = false;


  model = Random; // m
  numNodes = 100000; //n
  numEdges = 1000000; // e
  edgeProbability = 0.0001; //ep
  
  rmat.a = 0.45; // a
  rmat.b = 0.15; // b
  rmat.c = 0.15; // c
  rmat.d = 0.25; // d


  /**** Ma Clique Size **********/
  maxCliqueSize = 4;

  /*** Properties of interest *******/
  degree  = 10;// dg
  density = 0.0001; // dn
  clusteringCoeff = 0.1; // cc
  aed = 0.33; // aed
  degreesd = 0.2; //dsd
  err = 1; /*1%*/ // err
  maxWeight = 100; //mw
  minWeight = 1; // iw
}

bool GraphProperty::updateConfigs(const std::string configfile) {
  FILE *f;

  f = fopen(configfile.c_str(), "r");
  if(f == NULL) {
    printError(CONFIG_FILE_NOT_FOUND, 0, NULL);
    return false;
  }
  
  
  // we assume the first entry in the config file specifies the model
  // of graph generation.
  
  
  char fileType[100];
  
  fscanf(f,"%s", fileType);

  /*
    TODO allow multiple graph property input options
    right now we allow only allproperty option.
  */
  if(strcmp(fileType,"allproperty") == 0) {
    property = allproperty;
    while (!feof (f)) {
      char propType[50];
      double val;
      fscanf(f, "%s %lf", propType, &val);
      if(strcmp(propType,"degree") == 0) {
	degree = val;
      } else if(strcmp(propType,"density") == 0) {
	density = val;
      } else if(strcmp(propType,"cc") == 0) {
	clusteringCoeff = val;
      } else if(strcmp(propType,"aed") == 0) {
	aed = val;
      } else if(strcmp(propType,"dsd") == 0) {
	degreesd = val;
      } 
    }
  } else {

  }
  
}




