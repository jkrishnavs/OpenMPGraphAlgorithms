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
  edgeProbability = 0.00; //ep
  
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
  char modelType[100];
  
  int t = fscanf(f,"%s", fileType);
  
  t = fscanf(f,"%s", modelType);

  assert(t!=0);

  if(!strcmp(modelType, "random")) {
    model = Random;
  } else if(!strcmp(modelType, "ErdosRenyi")) {
    model = ErdosRenyi;
  } else if(!strcmp(modelType, "rmat")) {
    model = RMAT;
  } else if(!strcmp(modelType, "auto")) {
    model = Auto;
  } else{
    return false;
  }

  /*
    TODO allow multiple graph property input options
    right now we allow only allproperty option.
  */
  if(!strcmp(fileType,"allproperty")) {
    property = allproperty;
    while (!feof (f)) {
      char propType[50];
      double val;
      int t = fscanf(f, "%s %lf", propType, &val);
      assert(t!=0);
      if(!strcmp(propType,"degree")) {
	degree = val;
      } else if(!strcmp(propType,"density")) {
	density = val;
      } else if(!strcmp(propType,"cc")) {
	clusteringCoeff = val;
      } else if(!strcmp(propType,"aed")) {
	aed = val;
      } else if(!strcmp(propType,"dsd")) {
	degreesd = val;
      } else {
	fclose(f);
	return false;
      }
    }
  } else if(!strcmp(fileType, "degreeanddensity")) {
    // We use random generator
    property=degreeanddensity;
    while (!feof (f)) {
      char propType[50];
      double val;
      int t = fscanf(f, "%s %lf", propType, &val);
      assert(t!=0);
      if(!strcmp(propType,"degree")) {
	degree = val;
      } else if(!strcmp(propType,"density")) {
	density = val;
      } else {
	fclose(f);
	return false;
      }
      numNodes =  (node_t)  degree/density;
      numEdges = (edge_t) numNodes * degree;
    }
  } else{
    printf("Hello22");

    fclose(f);
    return false;
  }

  fclose(f);
  return true;
}




