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



void setdefaultValues() {
  selfloop = false;
  
}




int runalgo(int argc, char** argv) {
  
}
