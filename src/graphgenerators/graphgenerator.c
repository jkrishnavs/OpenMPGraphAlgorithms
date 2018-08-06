/*****
 * Original base algorithm from GTGraph. 
 ****/
#include<unistd.h>
#include"graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"




typedef enum GraphModel {
  Random,
  ErdosRenyi,
  RMAT,
  SSCA
} GraphModel;
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
  degree  = 10;// dg
  density = 0.0001; // dn
  clusteringCoeff = 0.1; // cc
  aed = 0.33; // aed
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

}




int runalgo(int argc, char** argv) {

  if(argc > 4) {
    const char* argList[3] = {" <configfile>" , "<outputfile>","<propfile>"};
    printError(INCORRECT_ARG_LIST, NO_OF_ARGS, argList);
    return -1;
  }

  
  
  
}
