/*****
 * Original base algorithm from GTGraph. 
 ****/
#include<unistd.h>
#include"graph.h"
#include "mainFunctions.h"
#include "print.h"
#include "powerperformacetracking.h"
#include "graphprop.h"

Graph* randomGenerator();
Graph* erdosRenyiGenerator();
Graph* rmatGenerator();
Graph* SSCAGenerator();
Graph* propertyControlledGraphGenerator();
Graph* callappropriategenerator();
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
  Graph * G = callappropriategenerator();
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

Graph* callappropriategenerator() {
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

Graph* randomGenerator() {
  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/
  LONG_T *startVertex, *endVertex;
  WEIGHT_T *weights;
  LONG_T i, j, u, v;
  WEIGHT_T w;
	LONG_T estNumEdges, numEdges, edgeNum;	
	int *stream1, *stream2;
	FILE* outfp;

	/*----------------------------------------------*/
	/*		initialize SPRNG 		*/
	/*----------------------------------------------*/

	stream1 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED1, SPRNG_DEFAULT);
	stream2 = init_sprng(SPRNG_CMRG, 0, 1, SPRNG_SEED2, SPRNG_DEFAULT);

	/*------------------------------------------------------*/
	/*		generate edges as per the		*/		 
	/*		graph model and user options	      	*/
	/*------------------------------------------------------*/

	if ((STORE_IN_MEMORY == 0) && (SORT_EDGELISTS == 0)) {
		fprintf(stderr, "Generating edges on the fly\n");
		outfp = fopen(OUTFILE, "w");
                fprintf(outfp, "c FILE			: %s\n", OUTFILE);
                fprintf(outfp, "c No. of vertices	: %ld\n", n);
		if (GRAPH_MODEL == 1) 
			fprintf(outfp, "c No. of edges		: %ld\n", m);
		else
			fprintf(outfp, "                                                 \n");
                fprintf(outfp, "c Max. weight	        : %ld\n", MAX_WEIGHT);
                fprintf(outfp, "c Min. weight          	: %ld\n", MIN_WEIGHT);
                fprintf(outfp, "c A directed arc from u to v of weight w\n");
                fprintf(outfp, "c is represented below as ' a  u  v  w '\n");
		fprintf(stderr, "Generating Edges ... ");
		
		/* Erdos-Renyi */ 	
		if (GRAPH_MODEL == 0) {
			
			/* Write the no. of edges later */
			fprintf(outfp, "                             \n");
			numEdges = 0;
			
			for (i=0; i<n; i++) {
				for (j=0; j<n; j++) {
					
					if ((i==j) && (SELF_LOOPS == 0))		
						continue;

					if (p > sprng(stream1)) {
						w = MIN_WEIGHT + (WEIGHT_T) \
						   (MAX_WEIGHT - MIN_WEIGHT) * sprng(stream2);
						/* Create edge */
						fprintf(outfp, "a %ld %ld %ld\n", i+1, j+1, w);
						numEdges++;
					}
							
				} 
			}

			m = numEdges;	
			fclose(outfp);
			fprintf(stderr, "done\n");
			updateEdgeVal(OUTFILE, numEdges);			

		} else {

                        fprintf(outfp, "p sp %ld %ld\n", n, m);
                        numEdges = 0;

                        for (i=0; i<m; i++) {

				u = (LONG_T) isprng(stream1) % n;
				v = (LONG_T) isprng(stream1) % n;
				if ((u == v) && (SELF_LOOPS == 0)) {
					i--;
					continue;
				}

				w = MIN_WEIGHT + (WEIGHT_T) (MAX_WEIGHT - MIN_WEIGHT) * sprng(stream2);

                                /* Create edge */
                                fprintf(outfp, "a %ld %ld %ld\n", u+1, v+1, w);
                               
                        }

			fclose(outfp);
                        fprintf(stderr, "done\n");
	
		}	

		free(stream1);
		free(stream2);
		return;

	}

	fprintf(stderr, "Generating edges ... ");

	if (GRAPH_MODEL == 0) {

		/* Estimate the no. of edges */
		if (SELF_LOOPS)
			estNumEdges = (LONG_T) (120 * n * n * p)/100;
		else
			estNumEdges = (LONG_T) (120 * n * (n-1) *  p)/100; 

		edgeNum = 0;
		numEdges = 0;

		startVertex = (LONG_T *) malloc(estNumEdges * sizeof(LONG_T)); 
		endVertex   = (LONG_T *) malloc(estNumEdges * sizeof(LONG_T));
		
		for (i=0; i<n; i++) {
			
			for (j=0; j<n; j++) {

				if ((i==j) && (SELF_LOOPS == 0))
                                	continue;

                                if (p > sprng(stream1)) {

					startVertex[edgeNum] = i; 
					endVertex[edgeNum]   = j;
					edgeNum++;

				}	
			}

		}

		numEdges = edgeNum;

	} else {

		startVertex = (LONG_T *) malloc(m * sizeof(LONG_T));
		endVertex   = (LONG_T *) malloc(m * sizeof(LONG_T));

		for (i=0; i<m; i++) {

			u = (LONG_T) isprng(stream1) % n;
               		v = (LONG_T) isprng(stream1) % n;
                	if ((u == v) && (SELF_LOOPS == 0)) {
				i--;
                		continue;
			}
			startVertex[i] = u;
			endVertex[i] = v;
		}		
		
		numEdges = m;	

	}

	fprintf(stderr, "done\n");

	free(stream1);

	/*----------------------------------------------*/
        /*              generate edge weights           */
        /*----------------------------------------------*/

        fprintf(stderr, "Generating edge weights ... ");

        weights = (WEIGHT_T *) malloc(numEdges*sizeof(WEIGHT_T));

        for (i=0; i<numEdges; i++) {
                weights[i] = MIN_WEIGHT + (WEIGHT_T) (MAX_WEIGHT - MIN_WEIGHT) \
                        * sprng(stream2);
        }

        fprintf(stderr, "done\n");
        free(stream2);

        /*-------------------------------------------------------*/
        /*              sort the edge lists with start           */
        /*              vertex as primary key                    */
        /*---------------------------------------------- --------*/

        if (SORT_EDGELISTS) {

                fprintf(stderr, "Sorting edge list by start vertex ... ");
		if (GRAPH_MODEL == 1) {
                	if (SORT_TYPE == 0) {
                        	/* Counting sort */
                        	countingSort(startVertex, endVertex, weights, numEdges);
                	} else {
                        	/* Heap sort */
                       		heapSort(startVertex, endVertex, weights, numEdges);
                	}
		}
                fprintf(stderr, "done\n");
        }

        g->start = startVertex;
        g->end = endVertex;
        g->w = weights;
        g->n = n;
        g->m = numEdges;
}

Graph* erdosRenyiGenerator() {

  /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/


}

Graph* rmatGenerator() {
    /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/

}


Graph* SSCAGenerator() {
    /**
   * NOTICE: The base algorithm follows GT Graph generators. For original
   * GTGraph code and  licence etc see GTgraph folder
   **/

  // TODO
}
